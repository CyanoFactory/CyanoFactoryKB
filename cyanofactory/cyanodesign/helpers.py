"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from PyNetMet2.enzyme import Enzyme
from PyNetMet2.metabolism import Metabolism
from bioparser.optgene import OptGeneParser

def apply_commandlist(model, commandlist):

    # Supported commands:
    # reaction, add, name, reaction, min_const, max_const
    # reaction, edit, name, reaction, min_const, max_const
    # (metabolite, add, name)
    # metabolite, edit, name, new_name, is_external

    """
    :type model: PyNetMet2.metabolism.Metabolism
    """
    for command in commandlist:
        if len(command) < 3:
            raise ValueError("Bad command " + str(command))

        if command[0] == "reaction":
            enzyme = model.get_reaction(command[2])

            if command[1] == "add":
                if len(command) < 7:
                    raise ValueError("Bad command " + str(command))

                if enzyme is not None:
                    raise ValueError("Reaction already in model: " + command[2])

                model.add_reaction(command[3])
                reaction = model.get_reaction(command[2])
                reaction.constraint = [command[4], command[5]]
                reaction.pathway = command[6]
                model.calcs()
            elif command[1] == "edit":
                if len(command) < 7:
                    raise ValueError("Bad command " + str(command))

                if enzyme is None:
                    raise ValueError("Reaction not in model: " + command[2])

                reaction = Enzyme(command[3])
                reaction.constraint = [command[4], command[5]]
                reaction.pathway = command[6]

                if reaction.name != enzyme.name and model.has_reaction(reaction.name):
                    raise ValueError("Reaction already in model: " + reaction.name)

                model.enzymes[model.dic_enzs[reaction.name]] = reaction
                model.calcs()
            elif command[1] == "delete":
                if enzyme is None:
                    raise ValueError("Reaction not in model: " + command[2])

                model.remove_reaction(enzyme.name)
            else:
                raise ValueError("Invalid operation " + command[1])

        elif command[0] == "metabolite":
            if command[1] in ["add", "edit", "delete"]:
                model.rename_metabolite(command[2], command[3])
                if command[4] and command[1] != "delete":
                    model.make_metabolite_external(command[3])
                else:
                    # These are deleted when non references them
                    model.make_metabolite_internal(command[3])
            else:
                raise ValueError("Invalid operation " + command[1])
        else:
            raise ValueError("Invalid command " + command[0])

    model.calcs()
    return model

def model_from_string(model_str):
    from bioparser.optgene import OptGeneParser
    from StringIO import StringIO
    import os
    import tempfile

    with tempfile.NamedTemporaryFile(delete=False) as fid:
        path = fid.name

        fid.write(model_str)

    try:
        org = Metabolism(path)
    except:
        return None
    finally:
        os.remove(path)

    return org

