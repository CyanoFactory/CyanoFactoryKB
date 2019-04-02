/*
Copyright (c) 2019 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
*/

class ElementBase {
    id: string;
    metaid: string;
    name: string;
    sbo_term: string;
    description: string;

    constructor() {
    }

    greet() {
        return "Hello, " + this.name;
    }

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
}

class MetabolicModel extends ElementBase {
    metabolites: Metabolite[];
    reactions: Reaction[];
    compartments: Compartment[];
    //objectives;
    //groups;

    static from_json(j: object) : MetabolicModel {
        let obj = new MetabolicModel();
        obj.read_attributes(j);
        return obj;
    }
}

class Compartment extends ElementBase {
    constant: boolean;
    units: string;

    read_attributes(attrs : object) {
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

    static from_json(j: object) : Compartment {
        let obj = new Compartment();
        obj.read_attributes(j);
        return obj;
    }
}

class Reaction extends ElementBase {
    lower_bound: number;
    upper_bound: number;
    lower_bound_name: string;
    upper_bound_name: string;
    reversible: boolean;
    fast: boolean;
    substrates: MetaboliteReference[];
    products: MetaboliteReference[];
    //gene_products;
    enabled: boolean;

    read_attributes(attrs : object) {
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
                } else if (k == "wedesign:lowerFlux") {
                    this.lower_bound = v;
                } else if (k == "wedesign:upperFlux") {
                    this.upper_bound = v;
                } else {
                    this.read_common_attribute(k, v);
                }
            }
        }
    }

    static from_json(j: object) : Reaction {
        let obj = new Reaction();
        obj.read_attributes(j);
        return obj;
    }
}

class MetaboliteReference extends ElementBase {
    constant: boolean;
    stoichiometry: number;

    read_attributes(attrs : object) {
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

    static from_json(j: object) : MetaboliteReference {
        let obj = new MetaboliteReference();
        obj.read_attributes(j);
        return obj;
    }
}

class Metabolite extends ElementBase {
    compartment: string;
    charge: number;
    formula: string;
    constant: boolean;
    boundary_condition: boolean;
    has_only_substance_units: boolean;

    read_attributes(attrs : object) {
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

    static from_json(j: object) : Metabolite {
        let obj = new Metabolite();
        obj.read_attributes(j);
        return obj;
    }
}
