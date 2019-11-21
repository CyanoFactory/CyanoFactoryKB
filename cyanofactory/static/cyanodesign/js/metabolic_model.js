/*
Copyright (c) 2019 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
*/
define(["require", "exports"], function (require, exports) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    var Internal;
    (function (Internal) {
        class ClassBuilder {
            constructor(typeObj) {
                this.typeObj = typeObj;
            }
            create() {
                return new this.typeObj();
            }
        }
        Internal.ClassBuilder = ClassBuilder;
        class LstOp {
            constructor(list, type) {
                this.lst = list;
                this.type = type;
            }
            get(key, value) {
                for (let item of this.lst) {
                    if (item[key] == value) {
                        return item;
                    }
                }
                return null;
            }
            checked_get(key, value) {
                let val = this.get(key, value);
                if (val == null) {
                    throw new Error(this.type.constructor.name + ": No object (" + key + ", " + value + ") found");
                }
                return val;
            }
            has(key, value) {
                for (const item of this.lst) {
                    if (item[key] == value) {
                        return true;
                    }
                }
                return false;
            }
            index(key, value) {
                let i = 0;
                for (const item of this.lst) {
                    if (item[key] == value) {
                        return i;
                    }
                    ++i;
                }
                return null;
            }
            checked_index(key, value) {
                let i = 0;
                for (const item of this.lst) {
                    if (item[key] == value) {
                        return i;
                    }
                    ++i;
                }
                throw new Error(this.type.constructor.name + ": No object (" + key + ", " + value + ") found");
            }
            add(item) {
                if (this.has("id", item.id)) {
                    throw new Error(this.type.constructor.name + ": " + item.id + " is already in the list");
                }
                this.lst.push(item);
            }
            remove(key, value) {
                let i = 0;
                let ret = false;
                while (i < this.lst.length) {
                    const item = this.lst[i];
                    if (item[key] == value) {
                        this.lst.splice(i, 1);
                        ret = true;
                    }
                    else {
                        ++i;
                    }
                }
                return ret;
            }
            create() {
                return new this.type();
            }
        }
        Internal.LstOp = LstOp;
    })(Internal = exports.Internal || (exports.Internal = {}));
    class ElementBase {
        constructor() {
            this.id = "";
            this.metaid = "";
            this.name = "";
            this.sbo_term = "";
            this.description = "";
        }
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
        get_name_or_id() {
            return this.name.length == 0 ? this.id : this.name;
        }
        fixup() { }
    }
    class Model extends ElementBase {
        constructor() {
            super(...arguments);
            this.metabolites = [];
            this.reactions = [];
            this.compartments = [];
            this.parameters = [];
            this.external_compartment = null;
            this.lower_bound_limit = 0;
            this.upper_bound_limit = 0;
        }
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
        readAttributes(attrs) {
            // no-op
        }
        fromJson(j) {
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
        }
        ;
        refreshConstraints() {
            let param_names = {};
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
        getDefaultCompartment() {
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
        refreshExternalCompartment() {
            this.external_compartment = null;
            return this.getExternalCompartment();
        }
        getExternalCompartment() {
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
        fba(glpk_worker, obj, maximize, create_exchange_reactions = false) {
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
            let transport_reacs = [];
            for (const met of this.metabolites) {
                subjectTo.push({
                    name: met.id,
                    vars: [],
                    bnds: { type: 5, ub: 0.0, lb: 0.0 }
                });
                if (create_exchange_reactions && met.isExternal(this)) {
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
                }
                else {
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
    exports.Model = Model;
    class Compartment extends ElementBase {
        constructor() {
            super(...arguments);
            this.constant = false;
            this.units = "";
        }
        updateId(new_id, model) {
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
    exports.Compartment = Compartment;
    class Reaction extends ElementBase {
        constructor() {
            super(...arguments);
            this.lower_bound = null;
            this.upper_bound = null;
            this.lower_bound_name = "";
            this.upper_bound_name = "";
            this.reversible = false;
            this.fast = false;
            this.substrates = [];
            this.products = [];
            //gene_products;
            this.enabled = true;
        }
        static fromBioOptString(bioopt_string, model) {
            function takeWhile(predicate, iterable, inclusive) {
                let arr = [];
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
            }
            let enzyme = new Reaction();
            let line = bioopt_string.split(/\s+/g);
            let res = takeWhile(function (x) {
                return x.indexOf(":") != x.length - 1;
            }, line, true);
            enzyme.name = res[0].join(" ").slice(0, -1).trim();
            enzyme.id = enzyme.name;
            line = res[1];
            if (line.length == 0) {
                throw { "enzyme": enzyme, "message": "Invalid reaction" };
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
                res = takeWhile(function (x) {
                    return x != "+" && x != "->" && x != "<->";
                }, line, false);
                // Parser is now on +, -> or <-> or EOL
                let mname = res[0].join(" ");
                line = res[1];
                // Test if there a multiple operators and reject them
                res = takeWhile(function (x) {
                    return x == "+" || x == "->" || x == "<->";
                }, line, false);
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
                        throw { "enzyme": enzyme, message: "More than one reaction arrow" };
                    }
                    arrow = line[0];
                    state = 1;
                }
                // Remove the operator
                line = line.slice(1);
            }
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
            enzyme.makeUnconstrained(model);
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
            this.read_list(j["listOfReactants"], this.substrates, new Internal.ClassBuilder(MetaboliteReference));
            this.read_list(j["listOfProducts"], this.products, new Internal.ClassBuilder(MetaboliteReference));
        }
        fixup() {
            const substrates = this.substrates;
            const products = this.products;
            this.substrates = [];
            this.products = [];
            this.read_list(substrates, this.substrates, new Internal.ClassBuilder(MetaboliteReference));
            this.read_list(products, this.products, new Internal.ClassBuilder(MetaboliteReference));
        }
        isConstrained(model) {
            if (this.reversible) {
                return this.lower_bound != model.lower_bound_limit || this.upper_bound != model.upper_bound_limit;
            }
            else {
                return this.lower_bound != 0 || this.upper_bound != model.upper_bound_limit;
            }
        }
        makeUnconstrained(model) {
            if (this.reversible) {
                this.lower_bound = model.lower_bound_limit;
            }
            else {
                this.lower_bound = 0.0;
            }
            this.upper_bound = model.upper_bound_limit;
        }
        updateId(new_id, model) {
            if (new_id == this.id) {
                return;
            }
            let idx = model.reaction.get("id", new_id);
            if (idx != null) {
                throw new Error("Reaction with ID " + new_id + " already exists");
            }
            this.id = new_id;
        }
        toHtml(model) {
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
        }
        ;
        reactionToString(model) {
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
        }
        ;
        toString(model) {
            let element = "";
            element = element.concat(this.get_name_or_id() + " : ");
            element = element.concat(this.reactionToString(model));
            return element;
        }
        constraintsToString(model) {
            if (!this.isConstrained(model)) {
                return "";
            }
            else {
                return "[" + this.lower_bound + ", " + this.upper_bound + "]";
            }
        }
        clearMetaboliteReference(model) {
            let affected_mets = [];
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
        }
        ;
        updateMetaboliteReference(model) {
            let affected_mets = [];
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
        }
        ;
        remove(model) {
            this.clearMetaboliteReference(model);
            let index = model.reactions.indexOf(this);
            if (index > -1) {
                model.reactions.splice(index, 1);
                return true;
            }
            return false;
        }
        ;
        getMetaboliteIds(model) {
            return this.substrates.concat(this.products).map(m => m.id);
        }
        ;
        getMetabolites(model) {
            return this.substrates.concat(this.products).map(m => model.metabolite.checked_get("id", m.id));
        }
        ;
    }
    exports.Reaction = Reaction;
    class MetaboliteReference extends ElementBase {
        constructor() {
            super(...arguments);
            this.constant = false;
            this.stoichiometry = 1;
        }
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
        getMetabolite(model, ignore_errors = false) {
            if (ignore_errors) {
                let met = model.metabolite.get("id", this.id);
                if (met == null) {
                    met = new Metabolite();
                    met.id = this.id;
                    met.name = this.id;
                }
                return met;
            }
            return model.metabolite.checked_get("id", this.id);
        }
        toHtml(model) {
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
        }
        ;
        toString() {
            return this.stoichiometry + " " + this.name;
        }
        ;
    }
    exports.MetaboliteReference = MetaboliteReference;
    class Metabolite extends ElementBase {
        constructor() {
            super(...arguments);
            this.compartment = "";
            this.charge = 0;
            this.formula = "";
            this.constant = false;
            this.boundary_condition = false;
            this.has_only_substance_units = false;
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
        updateId(new_id, model) {
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
        remove(model) {
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
        isUnused() {
            return this.consumed.length == 0 && this.produced.length == 0;
        }
        getReactions() {
            return this.consumed.concat(this.produced);
        }
        isExternal(model) {
            return this.compartment == model.getExternalCompartment().id;
        }
    }
    exports.Metabolite = Metabolite;
    class Parameter extends ElementBase {
        constructor() {
            super(...arguments);
            this.constant = false;
            this.units = "";
            this.value = 0;
        }
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
    exports.Parameter = Parameter;
});
