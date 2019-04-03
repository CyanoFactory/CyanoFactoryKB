/*
Copyright (c) 2019 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
*/
var MetabolicModel;
(function (MetabolicModel) {
    MetabolicModel.LOWER_BOUND_LIMIT = -100000.0;
    MetabolicModel.UPPER_BOUND_LIMIT = 100000.0;
    class ClassBuilder {
        constructor(typeObj) {
            this.typeObj = typeObj;
        }
        create() {
            return new this.typeObj();
        }
    }
    class LstOp {
        constructor(list) {
            this.lst = list;
        }
        get(key, value) {
            for (let item of this.lst) {
                if (item[key] == value) {
                    return item;
                }
            }
            return null;
        }
        has(key, value) {
            for (const item of this.lst) {
                if (item[key] == value) {
                    return true;
                }
            }
            return false;
        }
    }
    class Helper {
        static getById(id, lst) {
            for (const item of lst) {
                if (item.id == id) {
                    return item;
                }
            }
            return null;
        }
        static getByName(id, lst) {
            for (const item of lst) {
                if (item.name == id) {
                    return item;
                }
            }
            return null;
        }
    }
    class ElementBase {
        read_common_attribute(k, v) {
            if (k == "id") {
                this.id = v;
            }
            else if (k == "metaid") {
                this.metaid = v;
            }
            else if (k == "name") {
                this.name = v;
            }
            else if (k == "sboTerm") {
                this.sbo_term = v;
            }
        }
        read_list(src_list, dst_list, creator) {
            for (const item of src_list) {
                let obj = creator.create();
                obj.fromJson(item);
                dst_list.push(obj);
            }
        }
    }
    class Model extends ElementBase {
        constructor() {
            super(...arguments);
            this.metabolites = [];
            this.reactions = [];
            this.compartments = [];
            this.parameters = [];
        }
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
        readAttributes(attrs) {
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
            return this.metabolites.filter(function (value) {
                return value.compartment == "e";
            });
        }
        ;
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
    MetabolicModel.Model = Model;
    class Compartment extends ElementBase {
        readAttributes(attrs) {
            for (const k in attrs) {
                if (attrs.hasOwnProperty(k)) {
                    const v = attrs[k];
                    if (k == "constant") {
                        this.constant = v;
                    }
                    else if (k == "units") {
                        this.units = v;
                    }
                    else {
                        this.read_common_attribute(k, v);
                    }
                }
            }
        }
        fromJson(j) {
            this.readAttributes(j);
        }
    }
    MetabolicModel.Compartment = Compartment;
    class Reaction extends ElementBase {
        constructor() {
            super(...arguments);
            this.substrates = [];
            this.products = [];
            //gene_products;
            this.enabled = true;
        }
        static fromBioOptString(bioopt_string) {
            var takeWhile = function (predicate, iterable, inclusive) {
                var arr = [];
                iterable.every(function (x) {
                    if (predicate(x)) {
                        arr.push(x);
                    }
                    else {
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
            var res = takeWhile(function (x) { return x.indexOf(":") != x.length - 1; }, line, true);
            enzyme.name = res[0].join(" ").slice(0, -1).trim();
            enzyme.id = enzyme.name;
            line = res[1];
            if (line.length == 0) {
                throw { "enzyme": enzyme, "message": "Invalid reaction" };
            }
            var reactants = [];
            var arrow = "";
            var products = [];
            var state = 0;
            while (line.length > 0) {
                var value = 1;
                if (line[0].indexOf("/") != -1) {
                    var fl = line[0].split(/\//g);
                    value = +(fl[0]) / +(fl[1]);
                    if (isNaN(value)) {
                        value = 1;
                    }
                    else {
                        line = line.slice(1);
                    }
                }
                else {
                    value = +(line[0]);
                    if (isNaN(value)) {
                        value = 1;
                    }
                    else {
                        line = line.slice(1);
                    }
                }
                // Identifier...
                res = takeWhile(function (x) { return x != "+" && x != "->" && x != "<->"; }, line, false);
                // Parser is now on +, -> or <-> or EOL
                var mname = res[0].join(" ");
                line = res[1];
                // Test if there a multiple operators and reject them
                res = takeWhile(function (x) { return x == "+" && x == "->" && x == "<->"; }, line, false);
                // Parser is now on the next metabolite (number or name) or EOL
                if (res[0].length > 0) {
                    if (res[0].length > 1) {
                        throw { "enzyme": enzyme, message: "-> or <-> not allowed in names" };
                    }
                    // iname += ops[:-1]
                    // line.insert(0, ops[-1])
                    // Push operator back in line
                    // line.insert(0, ops[-1])
                }
                // Continue with old line, was just sanity check
                if (mname.length == 0) {
                    throw { "enzyme": enzyme, message: "Invalid reaction" };
                }
                let target;
                if (state == 0) {
                    target = enzyme.substrates;
                }
                else {
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
                        throw { "enzyme": enzyme, message: "More then one reaction arrow" };
                    }
                    arrow = line[0];
                    state = 1;
                }
                // Remove the operator
                line = line.slice(1);
            }
            ;
            if (arrow == "") {
                throw { "enzyme": enzyme, message: "No arrow found" };
            }
            if (enzyme.substrates.length == 0) {
                throw { "enzyme": enzyme, message: "No substrates found" };
            }
            if (enzyme.products.length == 0) {
                throw { "enzyme": enzyme, message: "No products found" };
            }
            enzyme.reversible = arrow == "<->";
            return enzyme;
        }
        ;
        static nameComparator(a, b) {
            if (a.name == null || b.name == null) {
                console.log(a);
            }
            return a.name.localeCompare(b.name);
        }
        ;
        readAttributes(attrs) {
            for (const k in attrs) {
                if (attrs.hasOwnProperty(k)) {
                    const v = attrs[k];
                    if (k == "fbc:lowerFluxBound") {
                        this.lower_bound_name = v;
                    }
                    else if (k == "fbc:upperFluxBound") {
                        this.upper_bound_name = v;
                    }
                    else if (k == "reversible") {
                        this.reversible = v;
                    }
                    else if (k == "fast") {
                        this.fast = v;
                    }
                    else if (k == "wedesign:enabled") {
                        this.enabled = v;
                    }
                    else {
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
        isConstrained() {
            return this.lower_bound != MetabolicModel.LOWER_BOUND_LIMIT && this.upper_bound != MetabolicModel.UPPER_BOUND_LIMIT;
        }
        makeUnconstrained() {
            this.lower_bound = MetabolicModel.LOWER_BOUND_LIMIT;
            this.upper_bound = MetabolicModel.UPPER_BOUND_LIMIT;
        }
        toHTML() {
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
        }
        ;
        reactionToString() {
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
        }
        ;
        toString() {
            var element = "";
            element = element.concat(this.name + " : ");
            element = element.concat(this.reactionToString());
            return element;
        }
        constraintsToString() {
            if (!this.isConstrained()) {
                return "";
            }
            else {
                return "[" + this.lower_bound + ", " + this.upper_bound + "]";
            }
        }
        clearMetaboliteReference() {
            var affected_mets = [];
            var that = this;
            // Clear refs
            MetabolicModel.model.metabolites.forEach(function (metabolite) {
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
        }
        ;
        updateMetaboliteReference() {
            var affected_mets = [];
            affected_mets = affected_mets.concat(this.clearMetaboliteReference());
            var that = this;
            // Add refs
            this.substrates.forEach(function (substrate) {
                var substrateObj = Helper.getById(substrate.id, MetabolicModel.model.metabolites);
                substrateObj.consumed.push(that);
                affected_mets.push(substrateObj);
                substrateObj.consumed = substrateObj.consumed.sort(Reaction.nameComparator);
            });
            this.products.forEach(function (product) {
                var productObj = Helper.getById(product.id, MetabolicModel.model.metabolites);
                productObj.produced.push(that);
                affected_mets.push(productObj);
                productObj.produced = productObj.produced.sort(Reaction.nameComparator);
            });
            return affected_mets;
        }
        ;
        remove() {
            this.clearMetaboliteReference();
            var index = MetabolicModel.model.reactions.indexOf(this);
            if (index > -1) {
                MetabolicModel.model.reactions.splice(index, 1);
                return true;
            }
            return false;
        }
        ;
        getMetabolites() {
            var nameProp = function (arg) {
                return arg.name;
            };
            var idGetter = function (arg) {
                Helper.getById(arg, MetabolicModel.model.metabolites);
            };
            return this.substrates.map(nameProp)
                .concat(this.products.map(nameProp))
                .map(idGetter);
        }
        ;
    }
    MetabolicModel.Reaction = Reaction;
    class MetaboliteReference extends ElementBase {
        readAttributes(attrs) {
            for (const k in attrs) {
                if (attrs.hasOwnProperty(k)) {
                    const v = attrs[k];
                    if (k == "constant") {
                        this.constant = v;
                    }
                    else if (k == "stoichiometry") {
                        this.stoichiometry = v;
                    }
                    else if (k == "species") {
                        this.id = v;
                    }
                    else {
                        this.read_common_attribute(k, v);
                    }
                }
            }
        }
        fromJson(j) {
            this.readAttributes(j);
        }
        toHtml() {
            var amount = this.stoichiometry + " ";
            var inst = Helper.getById(this.id, MetabolicModel.model.metabolites);
            let span = document.createElement("span");
            span.classList.add("cyano-metabolite");
            span.append(amount + inst.name);
            if (inst.isExternal()) {
                span.classList.add("cyano-external-metabolite");
            }
            return span;
        }
        ;
        toString() {
            return this.stoichiometry + " " + this.name;
        }
        ;
    }
    MetabolicModel.MetaboliteReference = MetaboliteReference;
    class Metabolite extends ElementBase {
        constructor() {
            super(...arguments);
            this.produced = [];
            this.consumed = [];
        }
        readAttributes(attrs) {
            for (const k in attrs) {
                if (attrs.hasOwnProperty(k)) {
                    const v = attrs[k];
                    if (k == "compartment") {
                        this.compartment = v;
                    }
                    else if (k == "fbc:charge") {
                        this.charge = v;
                    }
                    else if (k == "fbc:chemicalFormula") {
                        this.formula = v;
                    }
                    else if (k == "constant") {
                        this.constant = v;
                    }
                    else if (k == "boundaryCondition") {
                        this.boundary_condition = v;
                    }
                    else if (k == "hasOnlySubstanceUnits") {
                        this.has_only_substance_units = v;
                    }
                    else {
                        this.read_common_attribute(k, v);
                    }
                }
            }
        }
        fromJson(j) {
            this.readAttributes(j);
        }
        updateName(new_name) {
            let idx = Helper.getByName(new_name, MetabolicModel.model.metabolites);
            if (idx != null) {
                let that = this;
                this.consumed.concat(this.produced).forEach(function (enzyme) {
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
            this.consumed.forEach(function (enzyme) {
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
            var index = MetabolicModel.model.metabolites.indexOf(this);
            if (index > -1) {
                MetabolicModel.model.metabolites.splice(index, 1);
                return true;
            }
            return false;
        }
        isUnused() {
            return this.consumed.length == 0 && this.produced.length == 0;
        }
        getEnzymes() {
            return this.consumed.concat(this.produced);
        }
        isExternal() {
            return this.compartment == "e";
        }
    }
    MetabolicModel.Metabolite = Metabolite;
    class Parameter extends ElementBase {
        readAttributes(attrs) {
            for (const k in attrs) {
                if (attrs.hasOwnProperty(k)) {
                    const v = attrs[k];
                    if (k == "units") {
                        this.units = v;
                    }
                    else if (k == "constant") {
                        this.constant = v;
                    }
                    else if (k == "value") {
                        this.value = v;
                    }
                    else {
                        this.read_common_attribute(k, v);
                    }
                }
            }
        }
        fromJson(j) {
            this.readAttributes(j);
        }
    }
    MetabolicModel.Parameter = Parameter;
    MetabolicModel.model = new Model();
})(MetabolicModel || (MetabolicModel = {}));
