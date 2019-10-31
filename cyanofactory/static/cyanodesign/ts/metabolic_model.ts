/*
Copyright (c) 2019 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
*/

export namespace Internal {
    export class ClassBuilder {
        constructor(private typeObj: any) {

        }

        create() {
            return new this.typeObj();
        }
    }

    export class LstOp<T extends ElementBase> {
        readonly lst: T[];
        readonly type: any;

        constructor(list: T[], type: any) {
            this.lst = list;
            this.type = type;
        }

        get(key: string, value: any): T | null {
            for (let item of this.lst) {
                if (item[key] == value) {
                    return item;
                }
            }
            return null;
        }

        checked_get(key: string, value: any): T {
            let val = this.get(key, value);
            if (val == null) {
                throw new Error(this.type.constructor.name + ": No object (" + key + ", " + value + ") found");
            }
            return val;
        }

        has(key: string, value: any): boolean {
            for (const item of this.lst) {
                if (item[key] == value) {
                    return true;
                }
            }
            return false;
        }

        index(key: string, value: any): number | null {
            let i = 0;
            for (const item of this.lst) {
                if (item[key] == value) {
                    return i;
                }
                ++i;
            }
            return null;
        }

        checked_index(key: string, value: any): number {
            let i = 0;
            for (const item of this.lst) {
                if (item[key] == value) {
                    return i;
                }
                ++i;
            }
            throw new Error(this.type.constructor.name + ": No object (" + key + ", " + value + ") found");
        }

        add(item: T) {
            if (this.has("id", item)) {
                throw new Error(this.type.constructor.name + ": " + item.id + " is already in the list");
            }

            this.lst.push(item);
        }

        remove(key: string, value: any): boolean {
            let i = 0;
            let ret: boolean = false;
            while (i < this.lst.length) {
                const item: any = this.lst[i];
                if (item[key] == value) {
                    this.lst.splice(i, 1);
                    ret = true;
                } else {
                    ++i;
                }
            }
            return ret;
        }

        create(): T {
            return new this.type();
        }
    }
}

abstract class ElementBase {
    [key: string]: any;

    id: string = "";
    metaid: string = "";
    name: string = "";
    sbo_term: string = "";
    description: string = "";

    abstract fromJson(j: any): void;

    read_common_attribute(k: string, v: any) {
        if (k == "id") {
            this.id = v;
        } else if (k == "metaid") {
            this.metaid = v;
        } else if (k == "name") {
            this.name = v;
        } else if (k == "sboTerm") {
            this.sbo_term = v;
        }
    }

    read_list<T extends ElementBase>(src_list: any, dst_list: T[], creator: Internal.ClassBuilder) {
        for (const item of src_list) {
            let obj = creator.create();
            obj.fromJson(item);
            dst_list.push(obj);
        }
    }

    get_name_or_id(): string {
        return this.name.length == 0 ? this.id : this.name;
    }
}

export class Model extends ElementBase {
    metabolites: Metabolite[] = [];
    reactions: Reaction[] = [];
    compartments: Compartment[] = [];
    parameters: Parameter[] = [];

    external_compartment: Compartment = null;

    lower_bound_limit: number = 0;
    upper_bound_limit: number = 0;

    get metabolite() {
        return new Internal.LstOp(this.metabolites, Metabolite);
    }

    get reaction() {
        return new Internal.LstOp(this.reactions, Reaction);
    }

    get compartment() {
        return new Internal.LstOp(this.compartments, Compartment);
    }

    get parameter() {
        return new Internal.LstOp(this.parameters, Parameter);
    }

    //objectives;
    //groups;

    readAttributes(attrs: object) {
        // no-op
    }

    fromJson(j: any) {
        this.readAttributes(j);
        this.read_list(j["listOfSpecies"], this.metabolites, new Internal.ClassBuilder(Metabolite));
        this.read_list(j["listOfReactions"], this.reactions, new Internal.ClassBuilder(Reaction));
        this.read_list(j["listOfCompartments"], this.compartments, new Internal.ClassBuilder(Compartment));
        this.read_list(j["listOfParameters"], this.parameters, new Internal.ClassBuilder(Parameter));

        this.refreshConstraints();
    }

    getExternalMetabolites() {
        return this.metabolites.filter(function (value) {
            return value.compartment == "e";
        });
    };

    refreshConstraints() {
        let param_names: { [key: string]: Parameter } = {};
        for (let param of this.parameters) {
            param_names[param.id] = param;
        }

        for (let reaction of this.reactions) {
            reaction.lower_bound = param_names[reaction.lower_bound_name].value;
            this.lower_bound_limit = Math.min(this.lower_bound_limit, reaction.lower_bound);
            reaction.upper_bound = param_names[reaction.upper_bound_name].value;
            this.upper_bound_limit = Math.max(this.upper_bound_limit, reaction.upper_bound);
        }
    }

    getDefaultCompartment(): Compartment {
        if (this.compartments.length < 1) {
            return null;
        }

        if (this.compartments[0] == this.getExternalCompartment()) {
            if (this.compartments.length < 2) {
                return null;
            }
            return this.compartments[1];
        }

        return this.compartments[0];
    }

    refreshExternalCompartment(): Compartment {
        this.external_compartment = null;
        return this.getExternalCompartment();
    }

    getExternalCompartment(): Compartment {
        if (this.external_compartment != null) {
            return this.external_compartment;
        }

        // based on cobrapy
        const ext_list = ["e", "extracellular", "extraorganism", "out", "extracellular space",
          "extra organism", "extra cellular", "extra-organism", "external",
          "external medium"];

        for (const c of this.compartments) {
            for (const ext of ext_list) {
                if (c.id.toLowerCase() == ext) {
                    this.external_compartment = c;
                    return c;
                }
            }
        }

        for (const c of this.compartments) {
            if (c.id.toLowerCase().startsWith("e")) {
                this.external_compartment = c;
                return c;
            }
        }

        return null;
    }

    fba(glpk_worker: Worker, obj: Reaction, maximize: boolean) {
        let objective = {
            name: "z",
            direction: maximize ? 2 : 1,
            vars: [{
                name: obj.id,
                coef: 1.0
            }]
        };

        let subjectTo = [];
        let bounds = [];
        let transport_reacs: Reaction[] = [];

        for (const met of this.metabolites) {
            subjectTo.push({
                name: met.id,
                vars: [
                ],
                bnds: { type: 5, ub: 0.0, lb: 0.0 }
            });

            if (met.isExternal(this)) {
                let r = new Reaction();
                r.id = met.id + " <-> TRANSPORT";
                r.reversible = true;
                r.lower_bound = -10000.0;
                r.upper_bound = 10000.0;
                let mr = new MetaboliteReference();
                mr.id = met.id;
                mr.stoichiometry = 1.0;
                r.products.push(mr);
                transport_reacs.push(r);
            }
        }

        for (const reac of this.reactions) {
            for (const s of reac.substrates) {
                subjectTo[this.metabolite.index("id", s.id)]["vars"].push({
                    name: reac.id,
                    coef: -s.stoichiometry
                });
            }

            for (const p of reac.products) {
                subjectTo[this.metabolite.index("id", p.id)]["vars"].push({
                    name: reac.id,
                    coef: p.stoichiometry
                });
            }

            if (reac.enabled) {
                bounds.push({
                    name: reac.id,
                    type: 4,
                    ub: reac.upper_bound,
                    lb: reac.lower_bound
                });
            } else {
                bounds.push({
                    name: reac.id,
                    type: 5,
                    ub: 0.0,
                    lb: 0.0
                });
            }
        }

        for (const reac of transport_reacs) {
            for (const p of reac.products) {
                subjectTo[this.metabolite.index("id", p.id)]["vars"].push({
                    name: reac.id,
                    coef: p.stoichiometry
                });
            }

            bounds.push({
                name: reac.id,
                type: 4,
                ub: reac.upper_bound,
                lb: reac.lower_bound
            });
        }

        glpk_worker.postMessage({
            name: 'LP',
            objective: objective,
            subjectTo: subjectTo,
            bounds: bounds
        });
    }
}

export class Compartment extends ElementBase {
    constant: boolean = false;
    units: string = "";

    updateId(new_id: string, model: Model) {
        if (new_id == this.id) {
            return;
        }

        let idx = model.compartment.get("id", new_id);

        if (idx != null) {
            throw new Error("Compartment with ID " + new_id + " already exists");
        }

        for (let met of model.metabolites) {
            if (met.compartment == this.id) {
                met.compartment = new_id;
            }
        }

        this.id = new_id;
    }

    readAttributes(attrs: any) {
        for (const k in attrs) {
            if (attrs.hasOwnProperty(k)) {
                const v = attrs[k];
                if (k == "constant") {
                    this.constant = v;
                } else if (k == "units") {
                    this.units = v;
                } else {
                    this.read_common_attribute(k, v);
                }
            }
        }
    }

    fromJson(j: any) {
        this.readAttributes(j);
    }
}

export class Reaction extends ElementBase {
    lower_bound: number | null = null;
    upper_bound: number | null = null;
    lower_bound_name: string = "";
    upper_bound_name: string = "";
    reversible: boolean = false;
    fast: boolean = false;
    substrates: MetaboliteReference[] = [];
    products: MetaboliteReference[] = [];
    //gene_products;
    enabled: boolean = true;

    static fromBioOptString(bioopt_string: string): Reaction {
        function takeWhile(predicate: any, iterable: any, inclusive: any) {
            let arr: string[] = [];
            iterable.every(function (x: string) {
                if (predicate(x)) {
                    arr.push(x);
                } else {
                    if (inclusive) {
                        arr.push(x);
                    }
                    return false;
                }

                return true;
            });

            return [arr, iterable.slice(arr.length)];
        }

        let enzyme = new Reaction();

        let line = bioopt_string.split(/\s+/g);

        let res = takeWhile(function (x: string) {
            return x.indexOf(":") != x.length - 1
        }, line, true);

        enzyme.name = res[0].join(" ").slice(0, -1).trim();
        enzyme.id = enzyme.name;
        line = res[1];

        if (line.length == 0) {
            throw {"enzyme": enzyme, "message": "Invalid reaction"};
        }

        let reactants = [];
        let arrow = "";
        let products = [];
        let state = 0;

        while (line.length > 0) {
            let value = 1;

            if (line[0].indexOf("/") != -1) {
                let fl = line[0].split(/\//g);
                value = +(fl[0]) / +(fl[1]);
                if (isNaN(value)) {
                    value = 1;
                } else {
                    line = line.slice(1);
                }
            } else {
                value = +(line[0]);
                if (isNaN(value)) {
                    value = 1;
                } else {
                    line = line.slice(1);
                }
            }

            // Identifier...
            res = takeWhile(function (x: string) {
                return x != "+" && x != "->" && x != "<->"
            }, line, false);
            // Parser is now on +, -> or <-> or EOL
            let mname = res[0].join(" ");
            line = res[1];
            // Test if there a multiple operators and reject them
            res = takeWhile(function (x: string) {
                return x == "+" || x == "->" || x == "<->"
            }, line, false);

            // Parser is now on the next metabolite (number or name) or EOL
            if (res[0].length > 0) {
                if (res[0].length > 1) {
                    throw {"enzyme": enzyme, message: "-> or <-> not allowed in names"};
                }
                // iname += ops[:-1]
                // line.insert(0, ops[-1])
                // Push operator back in line
                // line.insert(0, ops[-1])
            }

            // Continue with old line, was just sanity check
            if (mname.length == 0) {
                throw {"enzyme": enzyme, message: "Invalid reaction"};
            }

            let target;
            if (state == 0) {
                target = enzyme.substrates;
            } else {
                target = enzyme.products;
            }

            let mr = new MetaboliteReference();
            mr.fromJson({
                "stoichiometry": value, "id": mname, "name": mname
            });
            target.push(mr);

            if (line.length == 0) {
                break;
            }

            if (line[0] == "->" || line[0] == "<->") {
                if (arrow != "") {
                    throw {"enzyme": enzyme, message: "More than one reaction arrow"};
                }

                arrow = line[0];
                state = 1;
            }

            // Remove the operator
            line = line.slice(1);
        }

        if (arrow == "") {
            throw {"enzyme": enzyme, message: "No arrow found"};
        }

        if (enzyme.substrates.length == 0) {
            throw {"enzyme": enzyme, message: "No substrates found"};
        }

        if (enzyme.products.length == 0) {
            throw {"enzyme": enzyme, message: "No products found"};
        }

        enzyme.reversible = arrow == "<->";
        enzyme.makeUnconstrained();

        return enzyme;
    };

    static nameComparator(a: Reaction, b: Reaction) {
        if (a.name == null || b.name == null) {
            console.log(a);
        }
        return a.name.localeCompare(b.name);
    };

    readAttributes(attrs: any) {
        for (const k in attrs) {
            if (attrs.hasOwnProperty(k)) {
                const v = attrs[k];
                if (k == "fbc:lowerFluxBound") {
                    this.lower_bound_name = v;
                } else if (k == "fbc:upperFluxBound") {
                    this.upper_bound_name = v;
                } else if (k == "reversible") {
                    this.reversible = v;
                } else if (k == "fast") {
                    this.fast = v;
                } else if (k == "wedesign:enabled") {
                    this.enabled = v;
                } else {
                    this.read_common_attribute(k, v);
                }
            }
        }
    }

    fromJson(j: any) {
        this.readAttributes(j);
        this.read_list(j["listOfReactants"], this.substrates, new Internal.ClassBuilder(MetaboliteReference));
        this.read_list(j["listOfProducts"], this.products, new Internal.ClassBuilder(MetaboliteReference));
    }

    isConstrained(model: Model): boolean {
        if (this.reversible) {
            return this.lower_bound != model.lower_bound_limit || this.upper_bound != model.upper_bound_limit;
        } else {
            return this.lower_bound != 0 || this.upper_bound != model.upper_bound_limit;
        }
    }

    makeUnconstrained(model: Model): void {
        this.lower_bound = model.lower_bound_limit;
        this.upper_bound = model.upper_bound_limit;
    }

    updateId(new_id: string, model: Model) {
        if (new_id == this.id) {
            return;
        }

        let idx = model.reaction.get("id", new_id);

        if (idx != null) {
            throw new Error("Reaction with ID " + new_id + " already exists");
        }

        this.id = new_id;
    }

    toHtml(model: Model): HTMLDivElement {
        let element = document.createElement("div");
        if (!this.enabled) {
            element.classList.add("cyano-reaction-disabled");
        }
        for (let i = 0; i < this.substrates.length; ++i) {
            element.appendChild(this.substrates[i].toHtml(model));
            if (i != this.substrates.length - 1) {
                element.appendChild(document.createTextNode(" + "));
            }
        }

        element.appendChild(document.createTextNode(this.reversible ? " ↔ " : " → "));

        for (let i = 0; i < this.products.length; ++i) {
            element.appendChild(this.products[i].toHtml(model));
            if (i != this.products.length - 1) {
                element.appendChild(document.createTextNode(" + "));
            }
        }

        return element;
    };

    reactionToString(model: Model): string {
        let element = "";

        for (let i = 0; i < this.substrates.length; ++i) {
            element = element.concat(this.substrates[i].stoichiometry + " ");

            element = element.concat(this.substrates[i].getMetabolite(model, true).get_name_or_id());
            if (i != this.substrates.length - 1) {
                element = element.concat(" + ");
            }
        }

        element = element.concat(this.reversible ? " <-> " : " -> ");

        for (let i = 0; i < this.products.length; ++i) {
            element = element.concat(this.products[i].stoichiometry + " ");

            element = element.concat(this.products[i].getMetabolite(model, true).get_name_or_id());
            if (i != this.products.length - 1) {
                element = element.concat(" + ");
            }
        }

        return element;
    };

    toString(model: Model): string {
        let element = "";
        element = element.concat(this.get_name_or_id() + " : ");
        element = element.concat(this.reactionToString(model));
        return element;
    }

    constraintsToString(model: Model): string {
        if (!this.isConstrained(model)) {
            return "";
        } else {
            return "[" + this.lower_bound + ", " + this.upper_bound + "]";
        }
    }

    clearMetaboliteReference(model: Model) {
        let affected_mets: Metabolite[] = [];

        let that = this;

        // Clear refs
        model.metabolites.forEach(function (metabolite) {
            let index = metabolite.consumed.indexOf(that);
            if (index > -1) {
                metabolite.consumed.splice(index, 1);
                affected_mets.push(metabolite);
                metabolite.consumed = metabolite.consumed.sort(Reaction.nameComparator);
            }

            index = metabolite.produced.indexOf(that);
            if (index > -1) {
                metabolite.produced.splice(index, 1);
                affected_mets.push(metabolite);
                metabolite.produced = metabolite.produced.sort(Reaction.nameComparator);
            }
        });

        return affected_mets;
    };

    updateMetaboliteReference(model: Model): Metabolite[] {
        let affected_mets: Metabolite[] = [];

        affected_mets = affected_mets.concat(this.clearMetaboliteReference(model));

        const that = this;

        // Add refs
        this.substrates.forEach(function (substrate) {
            let substrateObj = model.metabolite.checked_get("id", substrate.id);
            substrateObj.consumed.push(that);
            affected_mets.push(substrateObj);
            substrateObj.consumed = substrateObj.consumed.sort(Reaction.nameComparator);
        });

        this.products.forEach(function (product) {
            let productObj = model.metabolite.checked_get("id", product.id);
            productObj.produced.push(that);
            affected_mets.push(productObj);
            productObj.produced = productObj.produced.sort(Reaction.nameComparator);
        });

        return affected_mets;
    };

    remove(model: Model) {
        this.clearMetaboliteReference(model);

        let index = model.reactions.indexOf(this);
        if (index > -1) {
            model.reactions.splice(index, 1);
            return true;
        }
        return false;
    };

    getMetaboliteIds(model: Model): string[] {
        return this.substrates.concat(this.products).map(m => m.id);
    };

    getMetabolites(model: Model): Metabolite[] {
        return this.substrates.concat(this.products).map(
            m => model.metabolite.checked_get("id", m.id));
    };
}

export class MetaboliteReference extends ElementBase {
    constant: boolean = false;
    stoichiometry: number = 1;

    readAttributes(attrs: any) {
        for (const k in attrs) {
            if (attrs.hasOwnProperty(k)) {
                const v = attrs[k];
                if (k == "constant") {
                    this.constant = v;
                } else if (k == "stoichiometry") {
                    this.stoichiometry = v;
                } else if (k == "species") {
                    this.id = v;
                } else {
                    this.read_common_attribute(k, v);
                }
            }
        }
    }

    fromJson(j: any) {
        this.readAttributes(j);
    }

    getMetabolite(model: Model, ignore_errors: boolean = false): Metabolite {
        if (ignore_errors) {
            let met: Metabolite | null = model.metabolite.get("id", this.id);
            if (met == null) {
                met = new Metabolite();
                met.id = this.id;
                met.name = this.id;
            }
            return met;
        }

        return model.metabolite.checked_get("id", this.id);
    }

    toHtml(model: Model): HTMLSpanElement {
        let amount = this.stoichiometry + " ";
        let inst = this.getMetabolite(model);
        let span = document.createElement("span");
        span.classList.add("cyano-metabolite");
        span.dataset.id = this.id;
        span.append(amount + inst.name);

        if (inst.isExternal(model)) {
            span.classList.add("cyano-external-metabolite");
        }

        return span;
    };

    toString(): string {
        return this.stoichiometry + " " + this.name;
    };
}

export class Metabolite extends ElementBase {
    compartment: string = "";
    charge: number = 0;
    formula: string = "";
    constant: boolean = false;
    boundary_condition: boolean = false;
    has_only_substance_units: boolean = false;
    produced: Reaction[] = [];
    consumed: Reaction[] = [];

    readAttributes(attrs: any) {
        for (const k in attrs) {
            if (attrs.hasOwnProperty(k)) {
                const v = attrs[k];
                if (k == "compartment") {
                    this.compartment = v;
                } else if (k == "fbc:charge") {
                    this.charge = v;
                } else if (k == "fbc:chemicalFormula") {
                    this.formula = v;
                } else if (k == "constant") {
                    this.constant = v;
                } else if (k == "boundaryCondition") {
                    this.boundary_condition = v;
                } else if (k == "hasOnlySubstanceUnits") {
                    this.has_only_substance_units = v;
                } else {
                    this.read_common_attribute(k, v);
                }
            }
        }
    }

    fromJson(j: any) {
        this.readAttributes(j);
    }

    updateId(new_id: string, model: Model) {
        if (new_id == this.id) {
            return;
        }

        let idx = model.metabolite.get("id", new_id);

        if (idx != null) {
            throw new Error("Metabolite with ID " + new_id + " already exists");
        }

        let that = this;

        this.consumed.concat(this.produced).forEach(function (enzyme) {
            for (let enz of enzyme.substrates) {
                if (that.id == enz.id) {
                    enz.id = new_id;
                    break;
                }
            }

            for (let enz of enzyme.products) {
                if (that.id == enz.id) {
                    enz.id = new_id;
                    break;
                }
            }
        });

        this.id = new_id;
    }

    remove(model: Model) {
        let that = this;

        // Clear refs
        this.consumed.forEach(function (enzyme) {
            for (let i = 0; i < enzyme.substrates.length; ++i) {
                if (that.name == enzyme.substrates[i].name) {
                    enzyme.substrates.splice(i, 1);
                    break;
                }
            }

            for (let i = 0; i < enzyme.products.length; ++i) {
                if (that.name == enzyme.products[i].name) {
                    enzyme.products.splice(i, 1);
                    break;
                }
            }
        });

        let index = model.metabolites.indexOf(this);
        if (index > -1) {
            model.metabolites.splice(index, 1);
            return true;
        }
        return false;
    }

    isUnused(): boolean {
        return this.consumed.length == 0 && this.produced.length == 0;
    }

    getReactions(): Reaction[] {
        return this.consumed.concat(this.produced);
    }

    isExternal(model: Model): boolean {
        return this.compartment == model.getExternalCompartment().id;
    }
}

export class Parameter extends ElementBase {
    constant: boolean = false;
    units: string = "";
    value: number = 0;

    readAttributes(attrs: any) {
        for (const k in attrs) {
            if (attrs.hasOwnProperty(k)) {
                const v = attrs[k];
                if (k == "units") {
                    this.units = v;
                } else if (k == "constant") {
                    this.constant = v;
                } else if (k == "value") {
                    this.value = v;
                } else {
                    this.read_common_attribute(k, v);
                }
            }
        }
    }

    fromJson(j: any) {
        this.readAttributes(j);
    }
}

