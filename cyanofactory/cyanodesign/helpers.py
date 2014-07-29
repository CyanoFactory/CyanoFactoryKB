"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from PyNetMet.enzyme import Enzyme, EnzError
from PyNetMet.metabolism import Metabolism
from bioparser.optgene import OptGeneParser


def validate_bioopt(dic):
    """

    """

def delete_enzyme(model, name):
    pass

def rename_metabolite(model, met_from, met_to):
    OptGeneParser.Metabolite.instance(met_from).name = met_to

def make_metabolite_external(model, metabolite):
    OptGeneParser.Metabolite.instance(metabolite).external = True

def make_metabolite_internal(model, metabolite):
    OptGeneParser.Metabolite.instance(metabolite).external = False

def apply_commandlist(model, commandlist):
    # Supported commands:
    # reaction, add, name, reaction, min_const, max_const
    # reaction, edit, name, reaction, min_const, max_const
    # (metabolite, add, name)
    # metabolite, edit, name, new_name, is_external

    for command in commandlist:
        if len(command) < 3:
            raise ValueError("Bad command " + str(command))

        if command[0] == "reaction":
            enzyme = model.get_reaction(command[2])

            if command[1] == "add":
                if enzyme is not None:
                    raise ValueError("Reaction already in model: " + command[2])

                reaction = OptGeneParser.Reaction.from_string(command[3])
                model.add_reaction(reaction)
                reaction.constraint = [command[4], command[5]]

            elif command[1] == "edit":
                if enzyme is None:
                    raise ValueError("Reaction not in model: " + command[2])

                reaction = OptGeneParser.Reaction.from_string(command[3])
                reaction.constraint = [command[4], command[5]]

                if reaction.name != enzyme.name and model.has_reaction(reaction.name):
                    raise ValueError("Reaction already in model: " + reaction.name)

                enzyme.replace_with(reaction)
            elif command[1] == "delete":
                if enzyme is None:
                    raise ValueError("Reaction not in model: " + command[2])

                model.remove_reaction(enzyme)
            else:
                raise ValueError("Invalid operation " + command[1])

        elif command[0] == "metabolite":
            if command[1] == "edit":
                rename_metabolite(model, command[2], command[3])
                if command[4]:
                    make_metabolite_external(model, command[3])
                else:
                    make_metabolite_internal(model, command[3])
            else:
                raise ValueError("Invalid operation " + command[1])
        else:
            raise ValueError("Invalid command " + command[0])

    return model


def model_from_string(model_str):
    from bioparser.optgene import OptGeneParser
    from StringIO import StringIO

    org = OptGeneParser(StringIO(model_str))

    return org

