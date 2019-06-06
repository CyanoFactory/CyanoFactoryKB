"""
Copyright (c) 2019 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from metabolic_model import metabolic_model

def verify_object(obj, command):
    if "id" not in obj:
        raise ValueError("id not in command {}".format(command))

"""
Command format:
{
    "op": one of add, edit, delete
    "type": type of the SBML object modified
    "name": id of the SBML object modified
    "obj": changes applied to the object
"""

def apply_commandlist(model: metabolic_model.MetabolicModel, commandlist):
    """
    :type model: PyNetMet2.Metabolism
    """
    for command in commandlist:
        try:
            op = command["op"]
            typ = command["type"]
            idd = command["id"]
            obj = command["object"]
        except KeyError:
            raise ValueError("Bad command " + str(command))

        if typ not in ["reaction", "metabolite", "group"]:
            raise ValueError("Bad object type {} in command {}".format(typ, command))

        if op not in ["add", "edit", "delete"]:
            raise ValueError("Bad operation {} in command {}".format(op, command))

        if typ == "reaction":
            reaction = model.reaction.get(id=idd)  # type: metabolic_model.Reaction

            verify_object(obj, command)

            if op == "add":
                if reaction is not None:
                    raise ValueError("Reaction already in model: " + idd)
                if idd != obj["id"]:
                    raise ValueError("Reaction ID mismatch: " + idd)

                model.reactions.append(metabolic_model.Reaction.from_json(obj))
            elif op == "edit":
                if reaction is None:
                    raise ValueError("Reaction not in model: " + idd)

                if obj["id"] != reaction.id:
                    model.reaction.change_id(obj["id"], id=reaction.id)

                reaction.id = obj["id"]
                if "name" in obj:
                    reaction.name = obj["name"]

                if "substrates" in obj and "products" in obj:
                    reaction.substrates = list(
                        metabolic_model.MetaboliteReference.from_json(j) for j in obj["substrates"])
                    reaction.products = list(
                        metabolic_model.MetaboliteReference.from_json(j) for j in obj["products"])

                try:
                    reaction.reversible = obj["reversible"]
                except KeyError:
                    pass

                try:
                    reaction.lower_bound = obj["lower_bound"]
                    reaction.upper_bound = obj["upper_bound"]
                except KeyError:
                    pass

                # FIXME
                #try:
                #    reaction.pathway = obj["pathway"]
                #except KeyError:
                #    pass

                try:
                    reaction.enabled = obj["enabled"]
                except KeyError:
                    pass
            elif op == "delete":
                model.reaction.remove(id=reaction.id)
            else:
                raise ValueError("Invalid operation " + op)

        elif typ == "metabolite":
            if op in ["add", "edit"] and "id" not in obj:
                raise ValueError("Metabolite: Bad command " + str(command))

            if op == "add":
                if idd != obj["id"]:
                    raise ValueError("Metabolite ID mismatch: " + idd)

                model.metabolite.add(metabolic_model.Metabolite.from_json(obj))
            elif op == "edit":
                metabolite = model.metabolite.get(id=idd)
                if metabolite is None:
                    raise ValueError("Metabolite not in model: " + idd)

                try:
                    metabolite.compartment = obj["compartment"]
                except KeyError:
                    pass

                try:
                    metabolite.name = obj["name"]
                except KeyError:
                    pass

                if idd != obj["id"]:
                    model.metabolite.change_id(obj["id"], id=idd)

                    for reac in model.reactions:
                        for substrate in reac.substrates:
                            if substrate.id == idd:
                                substrate.id = obj["id"]
                        for product in reac.products:
                            if product.id == idd:
                                product.id = obj["id"]

            elif op == "delete":
                model.metabolite.remove(idd)
            else:
                raise ValueError("Invalid operation " + op)

        # FIXME
        #elif typ == "pathway":
        #    if "id" not in obj:
        #        raise ValueError("Bad command " + str(command))

        #    if op == "edit":
        #        for reac in model.enzymes:
        #            if reac.pathway == name:
        #                reac.pathway = obj["id"]

        #    else:
        #        raise ValueError("Invalid operation " + op)

        else:
            raise ValueError("Invalid command " + typ)

    return model

def compress_commandlist(commandlist):
    """
    :type model: json_model.JsonModel
    """
    last_command = None

    new_commandlist = []
    pending_commands = []

    # No error checking use after apply_commandlist

    for command in commandlist:
        op = command["op"]
        typ = command["type"]
        idd = command["id"]

        if len(pending_commands) == 0:
            pending_commands.append(command)
            last_command = command
            continue

        if typ == "metabolite" or typ == "reaction" or typ == "pathway":
            if last_command["type"] == typ:
                if last_command["object"].get("id", last_command["id"]) != idd or op == "add":
                    # Different to previous ones
                    new_commandlist += pending_commands
                    pending_commands = []
                    pending_commands.append(command)
                elif op == "edit":
                    if last_command["object"].get("id", last_command["id"]) == idd:
                        # Merge with previous
                        pending_commands[-1]["object"].update(command["object"])
                        # If previous was an add update the name
                        if pending_commands[-1]["op"] == "add":
                            pending_commands[-1]["id"] = command["object"]["id"]
                    else:
                        new_commandlist += pending_commands
                        pending_commands = []
                elif op == "delete":
                    # Remove all previous instances and take name of first
                    command["id"] = pending_commands[0]["id"]
                    # If first command was an "add" the item was created and just deleted
                    # can be removed completely
                    if pending_commands[0]["op"] != "add":
                        new_commandlist.append(command)

                    pending_commands = []
            else:
                new_commandlist += pending_commands
                pending_commands = []
                pending_commands.append(command)

        last_command = command

    new_commandlist += pending_commands
    return new_commandlist
