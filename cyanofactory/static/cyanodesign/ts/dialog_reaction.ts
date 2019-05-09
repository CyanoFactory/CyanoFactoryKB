import * as mm from "./metabolic_model";
import { DialogHelper, ElementWrapper } from "./dialog_helper";
import * as $ from "jquery";
import "datatables.net";
import {AppManager} from "./app";

declare var app : AppManager;

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
                            <label class="control-label" for="reaction-id">Identifier</label>
                            <input class="reaction-id form-control" type="text">
                        </div>
                        <div class="form-group">
                            <label class="control-label" for="reaction-name">Name</label>
                            <input class="reaction-name form-control" type="text">
                        </div>
                        <div class="form-group" style="display: none">
                            <label for="reaction-pathway">Pathway</label>
                            <select class="reaction-pathway" placeholder="Select or enter pathway"></select>
                        </div>
                        <div class="form-group">
                            <label for="reaction-substrates">Substrates</label>
                            <input class="reaction-substrates" type="text">
                        </div>
                        <div class="form-group">
                            <label for="reaction-products">Products</label>
                            <input class="reaction-products" type="text">
                        </div>
                        <div class="checkbox">
                            <input type="checkbox" name="enabled" class="enabled">
                            <label for="enabled">Enabled</label>
                        </div>
                        <div class="checkbox">
                            <input type="checkbox" name="reversible" class="reversible">
                            <label for="reversible">Reversible</label>
                        </div>
                        <div class="checkbox">
                            <input type="checkbox" name="constrained" class="constrained">
                            <label for="constrained">Constrained</label>
                        </div>
                        <div class="constraints-min-max form-inline">
                            <div class="form-group">
                                <label class="control-label" for="constrained-min">Min</label>
                                <input name="constrained_min" class="form-control constrained-min" type="text">
                            </div>
                            <div class="form-group">
                                <label class="control-label" for="constrained-max">Max</label>
                                <input name="constrained_max" class="form-control constrained-max" type="text">
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
    item: mm.Reaction = null;
    create: boolean = true;

    readonly id: ElementWrapper<string>;
    readonly name: ElementWrapper<string>;
    readonly substrates: ElementWrapper<string>;
    readonly products: ElementWrapper<string>;
    readonly enabled: ElementWrapper<boolean>;
    readonly reversible: ElementWrapper<boolean>;
    readonly constrained: ElementWrapper<boolean>;
    readonly constrained_min: ElementWrapper<string>;
    readonly constrained_max: ElementWrapper<string>;

    readonly substrate_refs: mm.MetaboliteReference[] = [];
    readonly product_refs: mm.MetaboliteReference[] = [];

    constructor(model: mm.Model) {
        document.body.appendChild(template.content.cloneNode(true));
        this.dialog_element = <HTMLElement>document.body.getElementsByClassName("dialog-reaction")[0]!;
        this.model = model;
        const self: Dialog = this;
        $(this.dialog_element).find(".btn-primary").click(function () {
            self.validate();
        });

        this.id = ElementWrapper.byClass<string>("reaction-id", this.dialog_element);
        this.name = ElementWrapper.byClass<string>("reaction-name", this.dialog_element);
        this.substrates = ElementWrapper.byClass<string>("reaction-substrates", this.dialog_element);
        this.products = ElementWrapper.byClass<string>("reaction-products", this.dialog_element);
        this.enabled = ElementWrapper.byClass<boolean>("enabled", this.dialog_element);
        this.reversible = ElementWrapper.byClass<boolean>("reversible", this.dialog_element);
        this.constrained = ElementWrapper.byClass<boolean>("constrained", this.dialog_element);
        this.constrained_min = ElementWrapper.byClass<string>("constrained-min", this.dialog_element);
        this.constrained_max = ElementWrapper.byClass<string>("constrained-max", this.dialog_element);

        for (const met of model.metabolites) {
            let sref = new mm.MetaboliteReference();
            let pref = new mm.MetaboliteReference();

            sref.id = met.id;
            sref.name = met.name;
            sref.stoichiometry = 1;

            pref.id = met.id;
            pref.name = met.name;
            pref.stoichiometry = 1;

            this.substrate_refs.push(sref);
            this.product_refs.push(pref);
        }
    }

    show(reaction: mm.Reaction) {
        if (reaction == null) {
            this.item = new mm.Reaction();
            this.create = true;
            this.item.makeUnconstrained();
        } else {
            this.item = reaction;
            this.create = false;
        }

        // Copy to dialog
        this.id.value = this.item.id;
        this.name.value = this.item.name;
        this.enabled.value = this.item.enabled;
        this.reversible.value = this.item.reversible;
        this.constrained.value = this.item.isConstrained();
        this.constrained_min.value = this.item.lower_bound.toString();
        this.constrained_max.value = this.item.upper_bound.toString();

        // Update dropdowns
        const obj_options: any = {
            valueField: 'id',
            labelField: 'name',
            searchField: ['name', 'id'],
            options: this.substrate_refs,
            persist: false,
            plugins: {
                'drag_drop' : {},
                'remove_button': {},
                'restore_on_backspace' : {
                    'text': function(item: mm.MetaboliteReference) {
                        return item.toString();
                    }
                }
            },
            render: {
                item: function(item: mm.MetaboliteReference, escape: any) {
                    return "<div>" + escape(item.toString()) + "</div>";
                },
                option: function(item: mm.MetaboliteReference, escape: any) {
                    return "<div>" + escape(item.get_name_or_id()) + "</div>";
                }
            },
            create: function (input: string) {
                let met = new mm.MetaboliteReference();
                met.id = input; // todo: sluggify
                met.name = input;
                met.stoichiometry = 1;

                const split: string[] = input.split(/\s+/, 1);
                if (split.length == 1) {
                    const num: number = DialogHelper.getFloat(split[0]);
                    if (Number.isFinite(num) && num > 0) {
                        met.stoichiometry = num;

                        let remain = input.substring(split[0].length + 1);
                        met.id = remain;
                        met.name = remain;
                    }
                }

                return met;
           }
        };

        $(this.substrates.element)["selectize"](obj_options);
        let substrates_selectize: any = (<any>$(this.substrates.element)[0]).selectize;
        substrates_selectize.clear();

        for (const metref of reaction.substrates) {
            substrates_selectize.addItem(metref.id);
        }

        $(this.products.element)["selectize"](obj_options);
        let products_selectize: any = (<any>$(this.products.element)[0]).selectize;
        products_selectize.clear();

        for (const metref of reaction.products) {
            products_selectize.addItem(metref.id);
        }

        // constraint div visibility
        const constraint_div: JQuery<HTMLElement> = $(this.dialog_element).find(".constraints-min-max");

        if (reaction.isConstrained()) {
            constraint_div.show();
        } else {
            constraint_div.hide();
        }

        $(this.dialog_element).find(".constrained:checkbox").change(function() {
            if ($(this).is(":checked")) {
                constraint_div.show();
            } else {
               constraint_div.hide();
            }
        });

        // clean up
        $(this.dialog_element).find("div").removeClass("has-error");
        $(this.dialog_element).find(".help-block").remove();

        // show
        $(this.dialog_element)["modal"]("show");
    }

    validate(): boolean {
        // clean up
        $(this.dialog_element).find("div").removeClass("has-error");
        $(this.dialog_element).find(".help-block").remove();

        // validate
        const constraints_enabled = $(this.dialog_element).find(".constrained").prop("checked");

        let valid = DialogHelper.checkId(this.id.element);
        valid = valid && DialogHelper.checkLength(this.name.element, "Name", 1, 255);
        valid = valid && DialogHelper.checkRegexp(this.name.element, /^[^:]+$/, "Name must not contain a colon (:)");

        let reac_with_id = this.model.reaction.get("id", this.id.value);

        if (reac_with_id != null) {
            valid = valid && DialogHelper.checkBool(this.id.element, reac_with_id == this.item, "Identifier already in use");
        }

        const cmin_float: number = DialogHelper.getFloat(this.constrained_min.value);
        const cmax_float: number = DialogHelper.getFloat(this.constrained_max.value);

        if (constraints_enabled) {
            const cmin_is_float: boolean = Number.isFinite(cmin_float);
            const cmax_is_float: boolean = Number.isFinite(cmax_float);
            valid = valid && DialogHelper.checkBool(this.constrained_min.element, cmin_is_float, "Min constraint must be numeric");
            valid = valid && DialogHelper.checkBool(this.constrained_min.element, cmax_is_float, "Max constraint must be numeric");
            if (cmin_is_float && cmax_is_float) {
                valid = valid && DialogHelper.checkBool(this.constrained_min.element, cmin_float <= cmax_float, "Min constraint must less or equal to Max constraint");
            }
        }

        if (!valid) {
            return false;
        }

        let reaction: mm.Reaction = this.item;
        if (this.create) {
            reaction.id = this.id.value;
            this.model.reactions.push(reaction);
        }

        let any_changed: boolean = this.create ||
            true || /* FIXME: substrates & products */
            reaction.id != this.id.value ||
            reaction.name != this.id.value ||
            reaction.enabled != this.enabled.value ||
            reaction.reversible != this.reversible.value ||
            reaction.isConstrained() != constraints_enabled ||
            reaction.lower_bound != cmin_float ||
            reaction.upper_bound != cmax_float;

        let old_id = reaction.id;
        reaction.updateId(this.id.value, this.model);

        reaction.name = this.name.value;
        reaction.enabled = this.enabled.value;
        reaction.reversible = this.reversible.value;

        /* constraints */
        if (!constraints_enabled) {
            reaction.makeUnconstrained();
        } else {
            reaction.lower_bound = cmin_float;
            reaction.upper_bound = cmax_float;
        }

        /* substrates & products */
        let substrates_selectize: any = (<any>$(this.substrates.element)[0]).selectize;
        let products_selectize: any = (<any>$(this.products.element)[0]).selectize;

        // remove all metabolites
        reaction.substrates = [];
        reaction.products = [];

        // add metabolites
        const met_adder = (mets: string[], target: mm.MetaboliteReference[]) => {
            for (const i of mets) {
                let met: mm.Metabolite | null = this.model.metabolite.get("id", i);
                if (met == null) {
                    // is a new metabolite
                    met = new mm.Metabolite();
                    met.id = i; // Todo: Sluggify
                    met.name = i;
                    app.metabolite_page.datatable.row.add(<any>met);
                    this.model.metabolites.push(met);
                }

                let metref = new mm.MetaboliteReference();
                metref.id = i;
                // FIXME!
                metref.stoichiometry = 1;

                target.push(metref);
            }
        };

        met_adder(substrates_selectize.items, reaction.substrates);
        met_adder(products_selectize.items, reaction.products);

        app.metabolite_page.datatable.sort();
        app.metabolite_page.datatable.draw();

        reaction.updateMetaboliteReference(this.model);

        if (any_changed) {
            app.command_list.push({
                "type": "reaction",
                "op": this.create ? "create" : "edit",
                "id": old_id,
                "object": {
                    "id": reaction.id,
                    "name": reaction.name,
                    "enabled": reaction.enabled,
                    "reversible": reaction.reversible,
                    "lower_bound": reaction.lower_bound,
                    "upper_bound": reaction.upper_bound,
                    "substrates": reaction.substrates.map((m: mm.MetaboliteReference) => m.id),
                    "products": reaction.products.map((m: mm.MetaboliteReference) => m.id)
                }
            });

            app.reaction_page.invalidate(reaction);
        }

        if (this.create) {
            app.reaction_page.datatable.row.add(reaction);
            app.reaction_page.datatable.sort();
            app.reaction_page.datatable.draw();
        }

        $(this.dialog_element)["modal"]("hide");

        return true;
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
