"use strict";

var Enzyme = (function() {
    function Enzyme(name) {
        this.name = name;
        this.reversible = false;
        this.substrates = [];
        this.products = [];
        this.stoichiometric_substrates = [];
        this.stoichiometric_products = [];
        this.constraints = [0, null];
        this.enabled = true;
        this.pathway = "GP";
    }

    Enzyme.prototype.isConstrained = function() {
        return this.constraints[1] !== null;
    };

    Enzyme.prototype.makeUnconstrained = function() {
        this.constraints[0] = this.reversible ? null : 0;
        this.constraints[1] = null;
    };

    Enzyme.prototype.toHTML = function() {
        var reac = function(substrates) {
            return function(value, index) {
                var amount = substrates[index] + " ";
                var name = $("<span>").addClass("cyano-metabolite").text(value.name);

                if (value.external) {
                    name.addClass("cyano-external-metabolite");
                }

                return amount + name.wrap("<p>").parent().html();
            }
        };

        var outer = $("<div></div>");
        var element = $("<div>");
        if (!this.enabled) {
            element.addClass("cyano-reaction-disabled");
        }
        element.append(this.substrates.map(reac(this.stoichiometric_substrates)).join(" + "));
        element.append(this.reversible ? " &harr; " : " &rarr; ");
        element.append(this.products.map(reac(this.stoichiometric_products)).join(" + "));

        element.append("</div>");
        element.appendTo(outer);

        return outer;
    };

    Enzyme.prototype.reactionToString = function() {
        var reac = function(substrates) {
            return function(value, index) {
                var elem = substrates[index] + " ";
                return elem + value.name;
            }
        };

        var element = "";
        element = element.concat(this.substrates.map(reac(this.stoichiometric_substrates)).join(" + "));
        element = element.concat(this.reversible ? " <-> " : " -> ");
        element = element.concat(this.products.map(reac(this.stoichiometric_products)).join(" + "));

        return element;
    };

    Enzyme.prototype.toString = function() {
        var element = "";
        element = element.concat(this.name + " : ");
        element = element.concat(this.reactionToString());
        return element;
    };

    Enzyme.prototype.constraintsToString = function() {
        if (!this.isConstrained()) {
            return "";
        } else {
            return "[" + this.constraints.join(", ") + "]";
        }
    };

    Enzyme.prototype.convertToFloat = function() {
        this.constraints[0] = parseFloat(this.constraints[0]);
        if (this.constraints[1] !== null) {
            this.constraints[1] = parseFloat(this.constraints[1]);
        }

        this.stoichiometric_substrates = this.stoichiometric_substrates.map(parseFloat);
        this.stoichiometric_products = this.stoichiometric_products.map(parseFloat);
    };

    Enzyme.prototype.clearMetaboliteReference = function() {
        var affected_mets = [];

        var that = this;

        // Clear refs
        Metabolite.metabolites.forEach(function(metabolite) {
            var index = metabolite.consumed.indexOf(that);
            if (index > -1) {
                metabolite.consumed.splice(index, 1);
                affected_mets.push(metabolite);
                metabolite.consumed = metabolite.consumed.sort(Metabolite.nameComparator);
            }

            index = metabolite.produced.indexOf(that);
            if (index > -1) {
                metabolite.produced.splice(index, 1);
                affected_mets.push(metabolite);
                metabolite.produced = metabolite.produced.sort(Metabolite.nameComparator);
            }
        });

        return affected_mets;
    };

    Enzyme.prototype.updateMetaboliteReference = function() {
        var affected_mets = [];

        affected_mets = affected_mets.concat(this.clearMetaboliteReference());

        var that = this;

        // Add refs
        this.substrates.forEach(function(substrate) {
            substrate.consumed.push(that);
            affected_mets.push(substrate);
            substrate.consumed = substrate.consumed.sort(Metabolite.nameComparator);
        });

        this.products.forEach(function(product) {
            product.produced.push(that);
            affected_mets.push(product);
            product.produced = product.produced.sort(Metabolite.nameComparator);
        });

        return affected_mets;
    };

    Enzyme.prototype.remove = function() {
        this.clearMetaboliteReference();

        var index = Enzyme.enzymes.indexOf(this);
        if (index > -1) {
            Enzyme.enzymes.splice(index, 1);
            return true;
        }
        return false;
    };

    Enzyme.prototype.getMetabolites = function() {
        return this.substrates.concat(this.products);
    };

    return Enzyme;
})();

Enzyme.enzymes = [];

Enzyme.indexByName = function(value) {
    for (var index = 0; index < Enzyme.enzymes.length; ++index) {
        if (Enzyme.enzymes[index].name == value) {
            return index;
        }
    }
    return -1;
};

Enzyme.fromReaction = function(reaction) {
    var enzyme = new Enzyme(reaction["name"]);

    enzyme.products = reaction.products.map(
        function(value) {
            return Metabolite.metabolites[Metabolite.indexByName(value)];
        }
    );

    enzyme.substrates = reaction.substrates.map(
        function(value) {
            return Metabolite.metabolites[Metabolite.indexByName(value)];
        }
    );

    enzyme.stoichiometric_substrates = reaction.stoichiometric[0];
    enzyme.stoichiometric_products = reaction.stoichiometric[1];
    enzyme.reversible = reaction.reversible;
    enzyme.constraints = reaction.constraints;
    enzyme.enabled = !reaction.disabled;
    enzyme.pathway = reaction.pathway;

    enzyme.updateMetaboliteReference();

    return enzyme;
};

Enzyme.getPathways = function() {
    var paths = Enzyme.enzymes.map(function(x) {
        return x.pathway;
    });

	var n = {};
    var r = [];
	for(var i = 0; i < paths.length; i++)
	{
		if (!n[paths[i]])
		{
			n[paths[i]] = true;
			r.push(paths[i]);
		}
	}
	return r.sort();
};

var Metabolite = (function() {
    function Metabolite(name, external) {
        this.external = external;
        this.name = name;
        this.consumed = [];
        this.produced = [];
    }

    Metabolite.prototype.updateName = function(new_name) {
        var idx = Metabolite.indexByName(new_name);

        if (idx === -1 ||
            Metabolite.metabolites[idx] == this) {
            this.name = new_name;
            Metabolite.metabolites[idx] = this;
            return true;
        }
        return false;
    };

    Metabolite.prototype.remove = function() {
        var that = this;

        // Clear refs
        this.consumed.forEach(function(enzyme) {
            var index = enzyme.substrates.indexOf(that);
            enzyme.substrates.splice(index, 1);
            enzyme.stoichiometric_substrates.splice(index, 1);

            var index = enzyme.products.indexOf(that);
            enzyme.products.splice(index, 1);
            enzyme.stoichiometric_products.splice(index, 1);
        });

        var index = Metabolite.metabolites.indexOf(this);
        if (index > -1) {
            Metabolite.metabolites.splice(index, 1);
            return true;
        }
        return false;
    };

    Metabolite.prototype.isUnused = function() {
        return this.consumed.length == 0 && this.produced == 0;
    };

    Metabolite.prototype.getEnzymes = function() {
        return this.consumed.concat(this.produced);
    };

    return Metabolite;
})();

Metabolite.metabolites = [];

Metabolite.indexByName = function(value) {
    for (var index = 0; index < Metabolite.metabolites.length; ++index) {
        if (Metabolite.metabolites[index].name == value) {
            return index;
        }
    }
    return -1;
};

Metabolite.fromReaction = function(reaction) {
    var metabolite_fn = function(value) {
        var idx = Metabolite.indexByName(value);
        if (idx == -1) {
            Metabolite.metabolites.push(new Metabolite(value, false));
        }
    };

    reaction["substrates"].forEach(metabolite_fn);
    reaction["products"].forEach(metabolite_fn);
};

Metabolite.createExternal = function(external) {
    var metabolite_fn = function(value) {
       var idx = Metabolite.indexByName(value);
        if (idx == -1) {
            Metabolite.metabolites.push(new Metabolite(value, true));
        } else {
            Metabolite.metabolites[idx].external = true;
        }
    };
    external.forEach(metabolite_fn);
};

Metabolite.getExternalMetabolites = function() {
    return Metabolite.metabolites.filter(function(value) {
        return value.external;
    });
};

Metabolite.nameComparator = function(a, b) {
    return a.name.localeCompare(b.name);
};
