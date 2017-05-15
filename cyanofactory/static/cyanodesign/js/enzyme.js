"use strict";

var Enzyme = (function() {
    function Enzyme(id, name) {
        this.id = id;
        this.name = name;
        this.substrates = [];
        this.products = [];
        this.constraints = [];
        this.reversible = false;
        this.pathway = "";
        this.disabled = false;
        this.favourite = false;
    }

    Enzyme.prototype.isConstrained = function() {
        return this.constraints.length > 0 && this.constraints[1] !== null;
    };

    Enzyme.prototype.makeUnconstrained = function() {
        this.constraints = []
    };

    Enzyme.prototype.toHTML = function() {
        var outer = $("<div></div>");
        var element = $("<div>");
        if (this.disabled) {
            element.addClass("cyano-reaction-disabled");
        }
        element.append(this.substrates.map(Enzyme.compoundToHTML).join(" + "));
        element.append(this.reversible ? " &harr; " : " &rarr; ");
        element.append(this.products.map(Enzyme.compoundToHTML).join(" + "));

        element.append("</div>");
        element.appendTo(outer);

        return outer;
    };

    Enzyme.prototype.reactionToString = function() {
        var element = "";
        element = element.concat(this.substrates.map(Enzyme.compoundToString).join(" + "));
        element = element.concat(this.reversible ? " <-> " : " -> ");
        element = element.concat(this.products.map(Enzyme.compoundToString).join(" + "));

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
        if (this.isConstrained()) {
            this.constraints[0] = parseFloat(this.constraints[0]);
            if (this.constraints[1] !== null) {
                this.constraints[1] = parseFloat(this.constraints[1]);
            }
        }

       var stoicFloat = function (value) {
            value.stoichiometry = parseFloat(value.stoichiometry);
            return value;
        };

        this.substrates = this.substrates.map(stoicFloat);
        this.products = this.products.map(stoicFloat);
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
            var substrateObj = Metabolite.instanceById(substrate.id);
            substrateObj.consumed.push(that);
            affected_mets.push(substrateObj);
            substrateObj.consumed = substrateObj.consumed.sort(Metabolite.nameComparator);
        });

        this.products.forEach(function(product) {
            var productObj = Metabolite.instanceById(product.id);
            productObj.produced.push(that);
            affected_mets.push(productObj);
            productObj.produced = productObj.produced.sort(Metabolite.nameComparator);
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
        var nameProp = function(arg) {
            return arg.name;
        };

        return this.substrates.map(nameProp).concat(this.products.map(nameProp)).map(Metabolite.instanceById);
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

Enzyme.compoundToHTML = function(compound) {
    var amount = compound.stoichiometry + " ";
    var inst = Metabolite.instanceById(compound.id);
    var name = $("<span>").addClass("cyano-metabolite").text(inst.name);

    if (inst.external) {
        name.addClass("cyano-external-metabolite");
    }

    return amount + name.wrap("<p>").parent().html();
};

Enzyme.compoundToString = function(compound) {
    return compound.stoichiometry  + " " + compound.name;
};

Enzyme.fromReaction = function(value) {
    var enzyme = new Enzyme(value.id, value.name);

    for (var key in value) {
        enzyme[key] = value[key];
    }

    enzyme.substrates.forEach(function(item) {
        item.name = Metabolite.instanceById(item.id).name;
    });
    enzyme.products.forEach(function(item) {
        item.name = Metabolite.instanceById(item.id).name;
    });

    return enzyme;
};

Enzyme.fromBioOptString = function(value) {
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

    var enzyme = new Enzyme();

    var line = value.split(/\s+/g);

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

        var target;
        if (state == 0) {
            target = enzyme.substrates;
        } else {
            target = enzyme.products;
        }

        target.push({"stoichiometry": value, "id": mname, "name": mname});

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

Enzyme.fromList = function(list) {
    for (var i = 0; i < list.length; ++i) {
        var enz = Enzyme.fromReaction(list[i]);
        Enzyme.enzymes.push(enz);
    }
};

var Metabolite = (function() {
    function Metabolite(id, name, external) {
        this.external = external;
        this.id = id;
        this.name = name;
        this.consumed = [];
        this.produced = [];
        this.chemical = "";
    }

    Metabolite.prototype.updateName = function(new_name) {
        var idx = Metabolite.indexByName(new_name);

        if (idx === -1 ||
            Metabolite.metabolites[idx] == this) {

            Metabolite.metabolites[idx] = this;

			var that = this;

			this.consumed.concat(this.produced).forEach(function(enzyme) {
				for (var i = 0; i < enzyme.substrates.length; ++i) {
					if (that.name == enzyme.substrates[i].name) {
						enzyme.substrates[i].name = new_name;
						break;
					}
				}

				for (var i = 0; i < enzyme.products.length; ++i) {
					if (that.name == enzyme.products[i].name) {
						enzyme.products[i].name = new_name
						break;
					}
				}
			});

			this.id = new_name;
			this.name = new_name;

            return true;
        }
        return false;
    };

    Metabolite.prototype.remove = function() {
        var that = this;

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

Metabolite.indexById = function(value) {
    for (var index = 0; index < Metabolite.metabolites.length; ++index) {
        if (Metabolite.metabolites[index].id == value) {
            return index;
        }
    }
    return -1;
};

Metabolite.indexByName = function(value) {
    for (var index = 0; index < Metabolite.metabolites.length; ++index) {
        if (Metabolite.metabolites[index].name == value) {
            return index;
        }
    }
    return -1;
};

Metabolite.instanceById = function(value) {
    for (var index = 0; index < Metabolite.metabolites.length; ++index) {
        if (Metabolite.metabolites[index].id == value) {
            return Metabolite.metabolites[index];
        }
    }
    return undefined;
};

Metabolite.instanceByName = function(value) {
    for (var index = 0; index < Metabolite.metabolites.length; ++index) {
        if (Metabolite.metabolites[index].name == value) {
            return Metabolite.metabolites[index];
        }
    }
    return undefined;
};

Metabolite.fromReaction = function(reaction) {
    var metabolite_fn = function(value) {
        var idx = Metabolite.indexById(value.name);
        if (idx == -1) {
            Metabolite.metabolites.push(new Metabolite(value.name, value.name, false));
        }
    };

    reaction["substrates"].forEach(metabolite_fn);
    reaction["products"].forEach(metabolite_fn);
};

Metabolite.fromList = function(list) {
    for (var i = 0; i < list.length; ++i) {
        Metabolite.metabolites.push(new Metabolite(list[i].id, list[i].name, list[i].external));
    }
};

Metabolite.createExternal = function(external) {
    var metabolite_fn = function(value) {
       var idx = Metabolite.indexByName(value);
        if (idx == -1) {
            Metabolite.metabolites.push(new Metabolite(value, value, true));
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
    if (a.name === undefined || b.name === undefined) {
        console.log(a);
    }
    return a.name.localeCompare(b.name);
};
