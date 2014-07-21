"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from PyNetMet.enzyme import Enzyme, EnzError
from PyNetMet.metabolism import Metabolism


def validate_bioopt(dic):
    """

    """


def get_enzyme_index(model, name):
    """
    :return: Index of enzyme name in model
    """
    for i, enzyme in enumerate(model.enzymes):
        if enzyme.name == name:
            return i

    return -1


def delete_enzyme(model, name):
    i = get_enzyme_index(model, name)



def rename_metabolite(model, met_from, met_to):
    for enzyme in model.enzymes:
        for i, metabolite in enumerate(enzyme.substrates):
            if metabolite == met_from:
                enzyme.substrates[i] = met_to

        for i, metabolite in enumerate(enzyme.products):
            if metabolite == met_from:
                enzyme.products[i] = met_to

    try:
        i = model.external.index(met_from)
        model.external[i] = met_to
    except ValueError:
        pass

    try:
        i = model.external_in.index(met_from)
        model.external_in[i] = met_to
    except ValueError:
        pass

    try:
        i = model.external_out.index(met_from)
        model.external_out[i] = met_to
    except ValueError:
        pass


def make_metabolite_external(model, metabolite):
    try:
        model.external.index(metabolite)
    except ValueError:
        model.external.append(metabolite)


def make_metabolite_internal(model, metabolite):
    try:
        model.external.remove(metabolite)
    except ValueError:
        pass

    try:
        model.external_in.remove(metabolite)
    except ValueError:
        pass

    try:
        model.external_out.remove(metabolite)
    except ValueError:
        pass


def apply_commandlist(model, commandlist):
    # After applying only calling dump is save because not all data is updated
    # properly

    # PyNetMet strips spaces everywhere, except in external... why?
    model.external = [ele.replace(" ","") for ele in model.external]

    # Supported commands:
    # reaction, add, name, reaction, min_const, max_const
    # reaction, edit, name, reaction, min_const, max_const
    # (metabolite, add, name)
    # metabolite, edit, name, new_name, is_external

    for command in commandlist:
        if len(command) < 3:
            raise ValueError("Bad command " + str(command))

        if command[0] == "reaction":
            enzyme_index = get_enzyme_index(model, command[2])

            if command[1] == "add":
                if enzyme_index != -1:
                    raise ValueError("Enzyme " + command[2] + " already in model")
                try:
                    Enzyme(command[3])
                except EnzError:
                    raise ValueError("Bad reaction " + command[3])
                model.add_reacs([command[3]])
                enzyme_index = get_enzyme_index(model, command[2])
                model.constr[enzyme_index] = (command[4], command[5])

            elif command[1] == "edit":
                try:
                    Enzyme(command[3])
                except EnzError:
                    raise ValueError("Bad reaction " + command[3])
                model.enzymes[enzyme_index] = Enzyme(command[3])
                model.constr[enzyme_index] = (command[4], command[5])
            elif command[1] == "delete":
                model.pop(enzyme_index)
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
    import os
    import tempfile
    
    with tempfile.NamedTemporaryFile(delete=False) as fid:
        path = fid.name

        fid.write(model_str)

    org = Metabolism(path)
    os.remove(path)

    return org

