import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

let template = document.createElement('template');
template.innerHTML = `
<div class="dialog-reaction modal fade" tabindex="-1" role="dialog" aria-labelledby="dialog-reaction-label"
     aria-hidden="true">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span
                        aria-hidden="true">&times;</span></button>
                <h4 class="modal-title" id="dialog-reaction-label">Edit Reaction</h4>
            </div>
            <div class="modal-body">
                <form>
                    <fieldset>
                        <div class="form-group">
                            <label class="control-label" for="enzyme-name">Name</label>
                            <input class="form-control" type="text" id="enzyme-name">
                        </div>
                        <div class="form-group" style="display: none">
                            <label for="reaction-pathway">Pathway</label>
                            <select id="reaction-pathway" placeholder="Select or enter pathway"></select>
                        </div>
                        <div class="form-group">
                            <label for="reaction-substrates">Substrates</label>
                            <input id="reaction-substrates" type="text">
                        </div>
                        <div class="form-group">
                            <label for="reaction-products">Products</label>
                            <input id="reaction-products" type="text">
                        </div>
                        <div class="checkbox">
                            <input type="checkbox" name="enabled" id="enabled">
                            <label for="enabled">Enabled</label>
                        </div>
                        <div class="checkbox">
                            <input type="checkbox" name="reversible" id="reversible">
                            <label for="reversible">Reversible</label>
                        </div>
                        <div class="checkbox">
                            <input type="checkbox" name="constrained" id="constrained">
                            <label for="constrained">Constrained</label>
                        </div>
                        <div id="constraints-min-max" class="form-inline">
                            <div class="form-group">
                                <label class="control-label" for="constrained-min">Min</label>
                                <input name="constrained_min" id="constrained-min" class="form-control" type="text">
                            </div>
                            <div class="form-group">
                                <label class="control-label" for="constrained-max">Max</label>
                                <input name="constrained_max" id="constrained-max" class="form-control" type="text">
                            </div>
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
    readonly model: mm.Model;
    readonly datatable: DataTables.Api;
    readonly dialog_element: HTMLElement;

    constructor(model: mm.Model) {
        document.body.appendChild(template.content.cloneNode(true));
        this.dialog_element = <HTMLElement>document.body.getElementsByClassName("dialog-reaction")[0]!;
        this.model = model;
    }

    show() {
        /*dialog_enzyme.data("object", enzyme);
        dialog_enzyme.data("create", create_new);
        dialog_enzyme.find("#enzyme-name").val(enzyme.id);
        dialog_enzyme.find("#reversible").prop("checked", enzyme.reversible);

        var substrates = dialog_enzyme.find("#reaction-substrates")[0].selectize;
        var products = dialog_enzyme.find("#reaction-products")[0].selectize;
        var enabled = dialog_enzyme.find("#enabled");
        var constrained = dialog_enzyme.find("#constrained");
        var constrained_min = dialog_enzyme.find("#constrained-min");
        var constrained_max = dialog_enzyme.find("#constrained-max");
        var pathway = dialog_enzyme.find("#reaction-pathway")[0].selectize;

        substrates.clear();
        products.clear();
        pathway.clear();

        substrates.clearOptions();
        products.clearOptions();
        pathway.clearOptions();

        model.metabolites.forEach(function(obj) {
            substrates.addOption({"id": obj.name, item: {obj: obj, amount: 1}});
            products.addOption({"id": obj.name, item: {obj: obj, amount: 1}});
        });
        substrates.refreshOptions();
        products.refreshOptions();

        var pathways = Enzyme.getPathways();
        for (var i = 0; i < pathways.length; ++i) {
            pathway.addOption({"text": pathways[i]});
        }
        pathway.refreshOptions();

        for (var i = 0; i < enzyme.substrates.length; ++i) {
            var substrate = enzyme.substrates[i];
            substrates.options[substrate.name].item.amount = substrate.stoichiometry;
            substrates.addItem(substrate.name);
        }
        for (i = 0; i < enzyme.products.length; ++i) {
            var product = enzyme.products[i];
            products.options[product.name].item.amount = product.stoichiometry;
            products.addItem(product.name);
        }

        substrates.refreshItems();
        products.refreshItems();
        pathway.addItem(enzyme.pathway);
        pathway.refreshItems();

        enabled.prop("checked", enzyme.enabled || create_new);

        constrained.prop("checked", enzyme.isConstrained());
        constrained_min.val(enzyme.constraints[0]);
        constrained_max.val(enzyme.constraints[1]);

        if (enzyme.isConstrained()) {
            constraints_div.show();
        } else {
            constraints_div.hide();
        }

        dialog_enzyme.find("div").removeClass("has-error");
        dialog_enzyme.find(".help-block").remove();*/

        $(this.dialog_element)["modal"]("show");
    }

    /*

        function saveEnzyme(enzyme) {
            var name = dialog_enzyme.find("#enzyme-name");
            var reversible = dialog_enzyme.find("#reversible");
            var enabled = dialog_enzyme.find("#enabled");
            var constrained = dialog_enzyme.find("#constrained");
            var constrained_min = dialog_enzyme.find("#constrained-min");
            var constrained_max = dialog_enzyme.find("#constrained-max");
            var substrates = dialog_enzyme.find("#reaction-substrates")[0].selectize;
            var products = dialog_enzyme.find("#reaction-products")[0].selectize;
            var pathway = dialog_enzyme.find("#reaction-pathway");

            var old_design_objective = design_objective;

            var visible_value = undefined;
            if (design_objective_visible_combobox[0].selectize.getValue().split('\x00').indexOf(enzyme.name) > -1) {
                visible_value = enzyme;
            }
            cyano_design_objective_select[0].selectize.removeOption(enzyme.name);
            cyano_design_design_objective_select[0].selectize.removeOption(enzyme.name);
            design_objective_visible_combobox[0].selectize.removeOption(enzyme.name);
            cyano_design_target_objective_select[0].selectize.removeOption(enzyme.name);

            if (enzyme.name != name.val()) {
                enzyme.id = name.val();
                enzyme.name = name.val();
            }

            cyano_design_objective_select[0].selectize.addOption(enzyme);
            cyano_design_objective_select[0].selectize.refreshOptions();
            cyano_design_design_objective_select[0].selectize.addOption(enzyme);
            cyano_design_design_objective_select[0].selectize.refreshOptions();
            design_objective_visible_combobox[0].selectize.addOption(enzyme);
            design_objective_visible_combobox[0].selectize.refreshOptions();
            cyano_design_target_objective_select[0].selectize.addOption(enzyme);
            cyano_design_target_objective_select[0].selectize.refreshOptions();

            if (visible_value !== undefined) {
                design_objective_visible_combobox[0].selectize.addItem(visible_value.name);
            }

            if (old_design_objective !== undefined) {
                cyano_design_objective_select[0].selectize.setValue(old_design_objective.name, false);
            }

            enzyme.reversible = reversible.prop("checked");
            if (constrained.prop("checked")) {
                enzyme.constraints[0] = constrained_min.val();
                enzyme.constraints[1] = constrained_max.val();
            } else {
                enzyme.makeUnconstrained();
            }

            enzyme.enabled = enabled.prop("checked");

            enzyme.substrates = [];
            enzyme.products = [];

            substrates.getValue().split('\x00').forEach(function(a) {
                enzyme.substrates.push({
                    "id": substrates.options[a].item.obj.id,
                    "name": substrates.options[a].item.obj.name,
                    "stoichiometry": substrates.options[a].item.amount
                });
            });

            products.getValue().split('\x00').forEach(function(a) {
                enzyme.products.push({
                    "id": products.options[a].item.obj.id,
                    "name": products.options[a].item.obj.name,
                    "stoichiometry": products.options[a].item.amount
                });
            });

            enzyme.pathway = pathway.val();

            enzyme.updateMetaboliteReference(model).forEach(function(m) {
                m.invalidate();
            });
            enzyme.invalidate();
        }

     */
}
