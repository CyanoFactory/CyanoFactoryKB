/*
Copyright (c) 2019 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
*/

namespace MetabolicModel {
    export const LOWER_BOUND_LIMIT = -100000.0;
    export const UPPER_BOUND_LIMIT = 100000.0;

    class ClassBuilder {
        constructor(private typeObj) {

        }
        create() {
            return new this.typeObj();
        }
    }

    class LstOp<T extends ElementBase> {
        private lst: T[];

        constructor(list: T[]) {
            this.lst = list;
        }

        get(key: string, value): T {
            for (let item of this.lst) {
                if (item[key] == value) {
                    return item;
                }
            }
            return null;
        }

        has(key: string, value): boolean {
            for (const item of this.lst) {
                if (item[key] == value) {
                    return true;
                }
            }
            return false;
        }
    }

    class Helper {
        static getById<T extends ElementBase>(id: string, lst: T[]): T {
            for (const item of lst) {
                if (item.id == id) {
                    return item;
                }
            }
            return null;
        }

        static getByName<T extends ElementBase>(id: string, lst: T[]): T {
            for (const item of lst) {
                if (item.name == id) {
                    return item;
                }
            }
            return null;
        }
    }

    abstract class ElementBase {
        id: string;
        metaid: string;
        name: string;
        sbo_term: string;
        description: string;

        abstract fromJson(j): void;

        read_common_attribute(k, v) {
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

        read_list<T extends ElementBase>(src_list, dst_list: T[], creator: ClassBuilder) {
            for (const item of src_list) {
                let obj = creator.create();
                obj.fromJson(item);
                dst_list.push(obj);
            }
        }
    }

    export class Model extends ElementBase {
        metabolites: Metabolite[] = [];
        reactions: Reaction[] = [];
        compartments: Compartment[] = [];
        parameters: Parameter[] = [];

        get metabolite() {
            return new LstOp(this.metabolites);
        }

        get reaction() {
            return new LstOp(this.reactions);
        }

        get compartment() {
            return new LstOp(this.compartments);
        }

        get parameter() {
            return new LstOp(this.parameters);
        }
        //objectives;
        //groups;

        readAttributes(attrs: object) {
            // no-op
        }

        fromJson(j) {
            this.readAttributes(j);
            this.read_list(j["listOfSpecies"], this.metabolites, new ClassBuilder(Metabolite));
            this.read_list(j["listOfReactions"], this.reactions, new ClassBuilder(Reaction));
            this.read_list(j["listOfCompartments"], this.compartments, new ClassBuilder(Compartment));
            this.read_list(j["listOfParameters"], this.parameters, new ClassBuilder(Parameter));

            this.refreshConstraints();
        }

        getExternalMetabolites() {
            return this.metabolites.filter(function(value) {
                return value.compartment == "e";
            });
        };

        refreshConstraints() {
            let param_names = {};
            for (let param of this.parameters) {
                param_names[param.id] = param;
            }

            for (let reaction of this.reactions) {
                reaction.lower_bound = param_names[reaction.lower_bound_name].value;
                reaction.upper_bound = param_names[reaction.upper_bound_name].value;
            }
        }
    }

    export class Compartment extends ElementBase {
        constant: boolean;
        units: string;

        readAttributes(attrs: object) {
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

        fromJson(j) {
            this.readAttributes(j);
        }
    }

    export class Reaction extends ElementBase {
        lower_bound: number;
        upper_bound: number;
        lower_bound_name: string;
        upper_bound_name: string;
        reversible: boolean;
        fast: boolean;
        substrates: MetaboliteReference[] = [];
        products: MetaboliteReference[] = [];
        //gene_products;
        enabled: boolean = true;

        static fromBioOptString(bioopt_string: string): Reaction {
            var takeWhile = function(predicate, iterable, inclusive) {
                var arr = [];
                iterable.every(function(x) {
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
            };

            var enzyme = new Reaction();

            var line = bioopt_string.split(/\s+/g);

            var res = takeWhile(function(x) { return x.indexOf(":") != x.length - 1 }, line, true);

            enzyme.name = res[0].join(" ").slice(0, -1).trim();
            enzyme.id = enzyme.name;
            line = res[1];

            if (line.length == 0) {
                throw {"enzyme": enzyme, "message": "Invalid reaction"};
            }

            var reactants = []
            var arrow = ""
            var products = []
            var state = 0

            while (line.length > 0) {
                var value = 1;

                if (line[0].indexOf("/") != -1) {
                    var fl = line[0].split(/\//g);
                    value = +(fl[0])/+(fl[1]);
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
                res = takeWhile(function(x) { return x != "+" && x != "->" && x != "<->" }, line, false);
                // Parser is now on +, -> or <-> or EOL
                var mname = res[0].join(" ");
                line = res[1];
                // Test if there a multiple operators and reject them
                res = takeWhile(function(x) { return x == "+" && x == "->" && x == "<->" }, line, false)

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

                target.push(new MetaboliteReference().fromJson({
                    "stoichiometry": value, "id": mname, "name": mname
                }));

                if (line.length == 0) {
                    break;
                }

                if (line[0] == "->" || line[0] == "<->") {
                    if (arrow != "") {
                        throw {"enzyme": enzyme, message: "More then one reaction arrow"};
                    }

                    arrow = line[0];
                    state = 1;
                }

                // Remove the operator
                line = line.slice(1);
            };

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

            return enzyme;
        };

        static nameComparator(a: Reaction, b: Reaction) {
            if (a.name == null || b.name == null) {
                console.log(a);
            }
            return a.name.localeCompare(b.name);
        };

        readAttributes(attrs: object) {
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

        fromJson(j) {
            this.readAttributes(j);
            this.read_list(j["listOfReactants"], this.substrates, new ClassBuilder(MetaboliteReference));
            this.read_list(j["listOfProducts"], this.products, new ClassBuilder(MetaboliteReference));
        }

        isConstrained(): boolean {
            return this.lower_bound != LOWER_BOUND_LIMIT && this.upper_bound != UPPER_BOUND_LIMIT;
        }

        makeUnconstrained(): void {
            this.lower_bound = LOWER_BOUND_LIMIT;
            this.upper_bound = UPPER_BOUND_LIMIT;
        }

        toHTML(): HTMLDivElement {
            let element = document.createElement("div");
            if (!this.enabled) {
                element.classList.add("cyano-reaction-disabled");
            }
            for (let i = 0; i < this.substrates.length; ++i) {
                element.appendChild(this.substrates[i].toHtml());
                if (i != this.substrates.length - 1) {
                    element.appendChild(document.createTextNode(" + "));
                }
            }

            element.appendChild(document.createTextNode(this.reversible ? " ↔ " : " → "));

            for (let i = 0; i < this.products.length; ++i) {
                element.appendChild(this.products[i].toHtml());
                if (i != this.products.length - 1) {
                    element.appendChild(document.createTextNode(" + "));
                }
            }

            return element;
        };

        reactionToString(): string {
            let element = "";

            for (let i = 0; i < this.substrates.length; ++i) {
                element.concat(this.substrates[i].toString());
                if (i != this.substrates.length - 1) {
                    element.concat(" + ");
                }
            }

            element = element.concat(this.reversible ? " <-> " : " -> ");

            for (let i = 0; i < this.products.length; ++i) {
                element.concat(this.products[i].toString());
                if (i != this.products.length - 1) {
                    element.concat(" + ");
                }
            }

            return element;
        };

        toString(): string {
            var element = "";
            element = element.concat(this.name + " : ");
            element = element.concat(this.reactionToString());
            return element;
        }

        constraintsToString(): string {
            if (!this.isConstrained()) {
                return "";
            } else {
                return "[" + this.lower_bound + ", " + this.upper_bound + "]";
            }
        }

        clearMetaboliteReference() {
            var affected_mets = [];

            var that = this;

            // Clear refs
            model.metabolites.forEach(function(metabolite) {
                var index = metabolite.consumed.indexOf(that);
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

        updateMetaboliteReference() {
            var affected_mets = [];

            affected_mets = affected_mets.concat(this.clearMetaboliteReference());

            var that = this;

            // Add refs
            this.substrates.forEach(function(substrate) {
                var substrateObj = Helper.getById(substrate.id, model.metabolites);
                substrateObj.consumed.push(that);
                affected_mets.push(substrateObj);
                substrateObj.consumed = substrateObj.consumed.sort(Reaction.nameComparator);
            });

            this.products.forEach(function(product) {
                var productObj = Helper.getById(product.id, model.metabolites);
                productObj.produced.push(that);
                affected_mets.push(productObj);
                productObj.produced = productObj.produced.sort(Reaction.nameComparator);
            });

            return affected_mets;
        };

        remove() {
            this.clearMetaboliteReference();

            var index = model.reactions.indexOf(this);
            if (index > -1) {
                model.reactions.splice(index, 1);
                return true;
            }
            return false;
        };

        getMetabolites() {
            var nameProp = function(arg) {
                return arg.name;
            };
            var idGetter = function(arg) {
                Helper.getById(arg, model.metabolites);
            };

            return this.substrates.map(nameProp)
                .concat(this.products.map(nameProp))
                .map(idGetter);
        };
    }

    export class MetaboliteReference extends ElementBase {
        constant: boolean;
        stoichiometry: number;

        readAttributes(attrs: object) {
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

        fromJson(j) {
            this.readAttributes(j);
        }

        toHtml(): HTMLSpanElement {
            var amount = this.stoichiometry + " ";
            var inst = Helper.getById(this.id, model.metabolites);
            let span = document.createElement("span")
            span.classList.add("cyano-metabolite");
            span.append(amount + inst.name);

            if (inst.isExternal()) {
                span.classList.add("cyano-external-metabolite");
            }

            return span;
        };

        toString(): string {
            return this.stoichiometry  + " " + this.name;
        };
    }

    export class Metabolite extends ElementBase {
        compartment: string;
        charge: number;
        formula: string;
        constant: boolean;
        boundary_condition: boolean;
        has_only_substance_units: boolean;
        produced: Reaction[] = [];
        consumed: Reaction[] = [];

        readAttributes(attrs: object) {
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

        fromJson(j) {
            this.readAttributes(j);
        }

        updateName(new_name) {
            let idx = Helper.getByName(new_name, model.metabolites);

            if (idx != null) {
                let that = this;

                this.consumed.concat(this.produced).forEach(function(enzyme) {
                    for (let i = 0; i < enzyme.substrates.length; ++i) {
                        if (that.name == enzyme.substrates[i].name) {
                            enzyme.substrates[i].name = new_name;
                            break;
                        }
                    }

                    for (let i = 0; i < enzyme.products.length; ++i) {
                        if (that.name == enzyme.products[i].name) {
                            enzyme.products[i].name = new_name;
                            break;
                        }
                    }
                });

                this.id = new_name;
                this.name = new_name;

                return true;
            }
            return false;
        }

        remove() {
            let that = this;

            // Clear refs
            this.consumed.forEach(function(enzyme) {
                for (var i = 0; i < enzyme.substrates.length; ++i) {
                    if (that.name == enzyme.substrates[i].name) {
                        enzyme.substrates.splice(i, 1);
                        break;
                    }
                }

                for (var i = 0; i < enzyme.products.length; ++i) {
                    if (that.name == enzyme.products[i].name) {
                        enzyme.products.splice(i, 1);
                        break;
                    }
                }
            });

            var index = model.metabolites.indexOf(this);
            if (index > -1) {
                model.metabolites.splice(index, 1);
                return true;
            }
            return false;
        }

        isUnused(): boolean {
            return this.consumed.length == 0 && this.produced.length == 0;
        }

        getEnzymes() {
            return this.consumed.concat(this.produced);
        }

        isExternal(): boolean {
            return this.compartment == "e";
        }
    }

    export class Parameter extends ElementBase {
        constant: boolean;
        units: string;
        value: number;

        readAttributes(attrs: object) {
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

        fromJson(j) {
            this.readAttributes(j);
        }
    }

    export let model = new Model();
}
