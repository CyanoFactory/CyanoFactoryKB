import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";
import {DialogHelper} from "./dialog_helper";

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

export class Dialog {
    app: app.AppManager;
    readonly datatable: DataTables.Api;
    readonly dialog_element: HTMLElement;
    readonly bulkdata_element: HTMLElement;
    readonly bulkdata_preview_element: HTMLElement;

    constructor(app: app.AppManager) {
        this.app = app;

        document.body.appendChild(template.content.cloneNode(true));
        this.dialog_element = <HTMLElement>document.body.getElementsByClassName("dialog-reaction-bulkadd")[0]!;
        this.bulkdata_element = <HTMLElement>document.body.getElementsByClassName("bulkdata")[0]!;
        this.bulkdata_preview_element = <HTMLElement>document.body.getElementsByClassName("bulkadd-preview")[0]!;

        const self: Dialog = this;
        $(this.dialog_element).find(".btn-primary").click(function () {
            self.validate();
        });

        // event handlers
        $(this.bulkdata_element).keyup(function() {
            const lines: string[] = (<string>$(this).val()).split(/\r*\n/);

            $(self.bulkdata_preview_element).empty();
            let buffered_line: string = "";

            for (const line of lines) {
                let e: string = "";
                try {
                    let last: string | null = line.split(/\s+/g).pop();
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
                } catch (err) {
                    e = err["enzyme"].toString() + " Error: " + err["message"];
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

    validate(): boolean {
        // clean up
        $(this.dialog_element).find("div").removeClass("has-error");
        $(this.dialog_element).find(".help-block").remove();

        const lines: string[] = (<string>$(this.bulkdata_element).val()).split(/\r*\n/);
        $(this.bulkdata_preview_element).empty();

        for (const line of lines) {
            //try {
                let mlen = this.app.model.metabolites.length;
                let reaction = mm.Reaction.fromBioOptString(line);



        }
            /*lines.forEach(function(line) {
                try {

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
        });*/

        $(this.dialog_element)["modal"]("hide");

        return true;
    }
}
