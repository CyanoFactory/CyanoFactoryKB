define(["require", "exports", "datatables.net"], function (require, exports) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    let template = document.createElement('template');
    template.innerHTML = `
<div class="dialog-reaction-bulkadd modal fade" tabindex="-1" role="dialog" aria-labelledby="dialog-reaction-bulkadd-label"
     aria-hidden="true">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span
                        aria-hidden="true">&times;</span></button>
                <h4 class="modal-title" id="dialog-reaction-bulkadd-label">Bulk add reactions</h4>
            </div>
            <div class="modal-body">
                <form>
                    <fieldset>
                        Enter reactions in BioOpt format for adding them to the model. Constraints can be put after the reactions.<br>
                        <div class="form-group">
                            <label class="control-label" for="enzyme-bulkdata">BioOpt data:</label>
                            <textarea class="form-control" id="enzyme-bulkdata" placeholder="reac1: a + 2 b -> 3.5 c [-1, 1]" rows="5"></textarea>
                        </div>
                        <div class="form-group" id="enzyme-bulkadd-preview">

                        </div>
                    </fieldset>
                </form>
            </div>
            <div class="modal-footer">
                <button type="button" class="btn btn-primary">Save changes</button>
                <button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
            </div>
        </div>
    </div>
</div>
`;
    class Dialog {
        constructor(model) {
            document.body.appendChild(template.content.cloneNode(true));
            this.dialog_element = document.body.getElementsByClassName("dialog-reaction-bulkadd")[0];
            this.model = model;
            /*
    
            $("#dialog-reaction-bulkadd").on("click", ".btn-primary", function(event) {
                var lines = $("#enzyme-bulkdata").val().split(/\r*\n/);
    
                $("#enzyme-bulkadd-preview").empty();
                lines.forEach(function(line) {
                    try {
                        var mlen = model.metabolites.length;
                        var reaction = Enzyme.fromBioOptString(line);
    
                        Metabolite.fromReaction(reaction);
                        if (model.metabolites.length > mlen) {
                            model.metabolites.slice(mlen).forEach(function (metabolite) {
                                command_list.push({
                                    "type": "metabolite",
                                    "op": "add",
                                    "id": metabolite.id,
                                    "object": {
                                        "id": metabolite.id,
                                        "name": metabolite.name,
                                        "external": metabolite.external
                                    }
                                });
    
                                datatable_metabolites.row.add(metabolite);
                            });
                        }
    
                        reaction.convertToFloat();
    
                        var command = {
                            "type": "reaction",
                            "id": reaction.id,
                            "object": jQuery.extend(true, {}, reaction)
                        };
    
                        if (Enzyme.indexByName(reaction.name) != -1) {
                            // Reaction exists
                            var orig_enzyme = model.reactions[Enzyme.indexByName(reaction.name)];
    
                            orig_enzyme.name = reaction.name;
                            orig_enzyme.substrates = reaction.substrates;
                            orig_enzyme.products = reaction.products;
                            orig_enzyme.reversible = reaction.reversible;
    
                            command["op"] = "edit";
                        } else {
    
                            model.reactions.push(reaction);
    
                            datatable_enzymes.row.add(reaction);
    
                            cyano_design_objective_select[0].selectize.addOption(reaction);
                            cyano_design_objective_select[0].selectize.refreshOptions();
                            cyano_design_design_objective_select[0].selectize.addOption(reaction);
                            cyano_design_design_objective_select[0].selectize.refreshOptions();
                            design_objective_visible_combobox[0].selectize.addOption(reaction);
                            design_objective_visible_combobox[0].selectize.refreshOptions();
                            cyano_design_target_objective_select[0].selectize.addOption(reaction);
                            cyano_design_target_objective_select[0].selectize.refreshOptions();
    
                            command["op"] = "add";
                        }
    
                        command_list.push(command);
    
                        reaction.updateMetaboliteReference(model).forEach(function(m) {
                            m.invalidate();
                        });
                        reaction.invalidate();
                    } catch (err) {
    
                    }
                });
    
                datatable_metabolites.sort();
                datatable_metabolites.draw();
    
                datatable_enzymes.sort();
                datatable_enzymes.draw();
    
                $("#dialog-reaction-bulkadd").modal("hide");
            });
    
             */
        }
    }
    exports.Dialog = Dialog;
});
