define(["require", "exports", "./metabolic_model", "jquery", "./metabolic_model", "datatables.net"], function (require, exports, mm, $, metabolic_model_1) {
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
                            <label class="control-label" for="bulkdata">BioOpt data:</label>
                            <textarea class="bulkdata form-control" placeholder="reac1: a + 2 b -> 3.5 c [-1, 1]" rows="5"></textarea>
                        </div>
                        <div class="bulkadd-preview form-group">

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
        constructor(app) {
            this.app = app;
            document.body.appendChild(template.content.cloneNode(true));
            this.dialog_element = document.body.getElementsByClassName("dialog-reaction-bulkadd")[0];
            this.bulkdata_element = document.body.getElementsByClassName("bulkdata")[0];
            this.bulkdata_preview_element = document.body.getElementsByClassName("bulkadd-preview")[0];
            const self = this;
            $(this.dialog_element).find(".btn-primary").click(function () {
                self.validate();
            });
            // event handlers
            $(this.bulkdata_element).keyup(function () {
                const lines = $(this).val().split(/\r*\n/);
                $(self.bulkdata_preview_element).empty();
                let buffered_line = "";
                for (const line of lines) {
                    let e = "";
                    try {
                        let last = line.split(/\s+/g).pop();
                        if (last == null) {
                            continue;
                        }
                        if (last == "<->" || last == "->" || last == "+") {
                            buffered_line += " " + line;
                            return;
                        }
                        buffered_line += " " + line;
                        buffered_line = buffered_line.trim();
                        if (buffered_line.length == 0) {
                            return;
                        }
                        e = mm.Reaction.fromBioOptString(buffered_line).toString(app.model);
                        buffered_line = "";
                    }
                    catch (err) {
                        e = err["enzyme"].toString(self.app.model) + " Error: " + err["message"];
                    }
                    if (e !== 'undefined') {
                        $("<div>").text(e).appendTo(".bulkadd-preview");
                    }
                }
            });
        }
        show() {
            // clean up
            $(this.dialog_element).find("div").removeClass("has-error");
            $(this.dialog_element).find(".help-block").remove();
            $(this.bulkdata_preview_element).empty();
            // show
            $(this.dialog_element)["modal"]("show");
        }
        validate() {
            // clean up
            $(this.dialog_element).find("div").removeClass("has-error");
            $(this.dialog_element).find(".help-block").remove();
            const lines = $(this.bulkdata_element).val().split(/\r*\n/);
            $(this.bulkdata_preview_element).empty();
            for (const line of lines) {
                try {
                    let mlen = this.app.model.metabolites.length;
                    let reaction = mm.Reaction.fromBioOptString(line);
                    for (const substrate of reaction.substrates) {
                        if (!this.app.model.metabolite.has("id", substrate.id)) {
                            let metabolite = new metabolic_model_1.Metabolite();
                            metabolite.id = substrate.id;
                            metabolite.name = substrate.name;
                            this.app.model.metabolites.push(metabolite);
                            this.app.command_list.push({
                                "type": "metabolite",
                                "op": "add",
                                "id": metabolite.id,
                                "object": {
                                    "id": metabolite.id,
                                    "name": metabolite.name,
                                    "external": false,
                                    "formula": metabolite.formula
                                }
                            });
                            this.app.metabolite_page.invalidate(metabolite);
                            this.app.metabolite_page.datatable.row.add(metabolite);
                        }
                    }
                    for (const product of reaction.products) {
                        if (!this.app.model.metabolite.has("id", product.id)) {
                            let metabolite = new metabolic_model_1.Metabolite();
                            metabolite.id = product.id;
                            metabolite.name = product.name;
                            this.app.model.metabolites.push(metabolite);
                            this.app.command_list.push({
                                "type": "metabolite",
                                "op": "add",
                                "id": metabolite.id,
                                "object": {
                                    "id": metabolite.id,
                                    "name": metabolite.name,
                                    "external": false,
                                    "formula": metabolite.formula
                                }
                            });
                            this.app.metabolite_page.invalidate(metabolite);
                            this.app.metabolite_page.datatable.row.add(metabolite);
                        }
                    }
                    this.app.metabolite_page.datatable.sort();
                    this.app.metabolite_page.datatable.draw();
                    reaction.makeUnconstrained();
                    let command = {
                        "type": "reaction",
                        "id": reaction.id,
                        "op": "add",
                        "object": {
                            "id": reaction.id,
                            "name": reaction.name,
                            "enabled": reaction.enabled,
                            "reversible": reaction.reversible,
                            "lower_bound": reaction.lower_bound,
                            "upper_bound": reaction.upper_bound,
                            "substrates": reaction.substrates,
                            "products": reaction.products
                        }
                    };
                    let orig_reac = this.app.model.reaction.get("id", reaction.id);
                    if (orig_reac != null) {
                        command["op"] = "edit";
                        orig_reac.substrates = reaction.substrates;
                        orig_reac.products = reaction.products;
                        orig_reac.reversible = reaction.reversible;
                        this.app.reaction_page.invalidate(orig_reac);
                        for (const met of orig_reac.updateMetaboliteReference(this.app.model)) {
                            this.app.metabolite_page.invalidate(met);
                        }
                    }
                    else {
                        this.app.model.reactions.push(reaction);
                        this.app.reaction_page.datatable.row.add(reaction);
                        for (const met of reaction.updateMetaboliteReference(this.app.model)) {
                            this.app.metabolite_page.invalidate(met);
                        }
                    }
                    this.app.command_list.push(command);
                    this.app.metabolite_page.datatable.draw();
                    this.app.reaction_page.datatable.sort();
                    this.app.reaction_page.datatable.draw();
                }
                catch (err) {
                }
            }
            this.app.settings_page.refresh();
            $(this.dialog_element)["modal"]("hide");
            return true;
        }
    }
    exports.Dialog = Dialog;
});
