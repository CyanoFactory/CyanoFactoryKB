define(["require", "exports", "./metabolic_model", "./dialog_helper", "jquery", "datatables.net", "selectize"], function (require, exports, mm, dialog_helper_1, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
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
    class Dialog {
        constructor(app) {
            this.item = null;
            this.create = true;
            this.substrate_refs = [];
            this.product_refs = [];
            this.app = app;
            document.body.appendChild(template.content.cloneNode(true));
            this.dialog_element = document.body.getElementsByClassName("dialog-reaction")[0];
            const self = this;
            $(this.dialog_element).find(".btn-primary").click(function () {
                self.validate();
            });
            this.id = dialog_helper_1.ElementWrapper.byClass("reaction-id", this.dialog_element);
            this.name = dialog_helper_1.ElementWrapper.byClass("reaction-name", this.dialog_element);
            this.substrates = dialog_helper_1.ElementWrapper.byClass("reaction-substrates", this.dialog_element);
            this.products = dialog_helper_1.ElementWrapper.byClass("reaction-products", this.dialog_element);
            this.enabled = dialog_helper_1.ElementWrapper.byClass("enabled", this.dialog_element);
            this.reversible = dialog_helper_1.ElementWrapper.byClass("reversible", this.dialog_element);
            this.constrained = dialog_helper_1.ElementWrapper.byClass("constrained", this.dialog_element);
            this.constrained_min = dialog_helper_1.ElementWrapper.byClass("constrained-min", this.dialog_element);
            this.constrained_max = dialog_helper_1.ElementWrapper.byClass("constrained-max", this.dialog_element);
            for (const met of this.app.model.metabolites) {
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
        show(reaction) {
            if (reaction == null) {
                this.item = new mm.Reaction();
                reaction = this.item;
                this.create = true;
                this.item.makeUnconstrained();
            }
            else {
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
            const obj_options = {
                valueField: 'id',
                labelField: 'name',
                searchField: ['name', 'id'],
                options: this.substrate_refs,
                persist: false,
                plugins: {
                    'drag_drop': {},
                    'remove_button': {},
                    'restore_on_backspace': {
                        'text': function (item) {
                            return item.toString();
                        }
                    }
                },
                render: {
                    item: function (item, escape) {
                        return "<div>" + escape(item.toString()) + "</div>";
                    },
                    option: function (item, escape) {
                        return "<div>" + escape(item.get_name_or_id()) + "</div>";
                    }
                },
                create: function (input) {
                    let met = new mm.MetaboliteReference();
                    met.id = input; // todo: sluggify
                    met.name = input;
                    met.stoichiometry = 1;
                    const split = input.split(/\s+/, 1);
                    if (split.length == 1) {
                        const num = dialog_helper_1.DialogHelper.getFloat(split[0]);
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
            $(this.substrates.element).selectize(obj_options);
            let substrates_selectize = this.substrates.element.selectize;
            substrates_selectize.clear();
            for (const metref of reaction.substrates) {
                substrates_selectize.addItem(metref.id);
            }
            $(this.products.element)["selectize"](obj_options);
            let products_selectize = this.products.element.selectize;
            products_selectize.clear();
            for (const metref of reaction.products) {
                products_selectize.addItem(metref.id);
            }
            // constraint div visibility
            const constraint_div = $(this.dialog_element).find(".constraints-min-max");
            if (reaction.isConstrained()) {
                constraint_div.show();
            }
            else {
                constraint_div.hide();
            }
            $(this.dialog_element).find(".constrained:checkbox").change(function () {
                if ($(this).is(":checked")) {
                    constraint_div.show();
                }
                else {
                    constraint_div.hide();
                }
            });
            // clean up
            $(this.dialog_element).find("div").removeClass("has-error");
            $(this.dialog_element).find(".help-block").remove();
            // show
            $(this.dialog_element)["modal"]("show");
        }
        validate() {
            // clean up
            $(this.dialog_element).find("div").removeClass("has-error");
            $(this.dialog_element).find(".help-block").remove();
            // validate
            const constraints_enabled = $(this.dialog_element).find(".constrained").prop("checked");
            let valid = dialog_helper_1.DialogHelper.checkId(this.id.element);
            valid = valid && dialog_helper_1.DialogHelper.checkLength(this.name.element, "Name", 1, 255);
            valid = valid && dialog_helper_1.DialogHelper.checkRegexp(this.name.element, /^[^:]+$/, "Name must not contain a colon (:)");
            let reac_with_id = this.app.model.reaction.get("id", this.id.value);
            if (reac_with_id != null) {
                valid = valid && dialog_helper_1.DialogHelper.checkBool(this.id.element, reac_with_id == this.item, "Identifier already in use");
            }
            const cmin_float = dialog_helper_1.DialogHelper.getFloat(this.constrained_min.value);
            const cmax_float = dialog_helper_1.DialogHelper.getFloat(this.constrained_max.value);
            if (constraints_enabled) {
                const cmin_is_float = Number.isFinite(cmin_float);
                const cmax_is_float = Number.isFinite(cmax_float);
                valid = valid && dialog_helper_1.DialogHelper.checkBool(this.constrained_min.element, cmin_is_float, "Min constraint must be numeric");
                valid = valid && dialog_helper_1.DialogHelper.checkBool(this.constrained_min.element, cmax_is_float, "Max constraint must be numeric");
                if (cmin_is_float && cmax_is_float) {
                    valid = valid && dialog_helper_1.DialogHelper.checkBool(this.constrained_min.element, cmin_float <= cmax_float, "Min constraint must less or equal to Max constraint");
                }
            }
            if (!valid) {
                return false;
            }
            let reaction = this.item;
            if (this.create) {
                reaction.id = this.id.value;
                this.app.model.reactions.push(reaction);
            }
            let any_changed = this.create ||
                true || /* FIXME: substrates & products */
                reaction.id != this.id.value ||
                reaction.name != this.id.value ||
                reaction.enabled != this.enabled.value ||
                reaction.reversible != this.reversible.value ||
                reaction.isConstrained() != constraints_enabled ||
                reaction.lower_bound != cmin_float ||
                reaction.upper_bound != cmax_float;
            let old_id = reaction.id;
            reaction.updateId(this.id.value, this.app.model);
            reaction.name = this.name.value;
            reaction.enabled = this.enabled.value;
            reaction.reversible = this.reversible.value;
            /* constraints */
            if (!constraints_enabled) {
                reaction.makeUnconstrained();
            }
            else {
                reaction.lower_bound = cmin_float;
                reaction.upper_bound = cmax_float;
            }
            /* substrates & products */
            let substrates_selectize = this.substrates.element.selectize;
            let products_selectize = this.products.element.selectize;
            // remove all metabolites
            reaction.substrates = [];
            reaction.products = [];
            // add metabolites
            const met_adder = (mets, target) => {
                for (const i of mets) {
                    let met = this.app.model.metabolite.get("id", i);
                    if (met == null) {
                        // is a new metabolite
                        met = new mm.Metabolite();
                        met.id = i; // Todo: Sluggify
                        met.name = i;
                        this.app.metabolite_page.datatable.row.add(met);
                        this.app.model.metabolites.push(met);
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
            this.app.metabolite_page.datatable.sort();
            this.app.metabolite_page.datatable.draw();
            reaction.updateMetaboliteReference(this.app.model);
            if (any_changed) {
                this.app.command_list.push({
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
                        "substrates": reaction.substrates.map((m) => m.id),
                        "products": reaction.products.map((m) => m.id)
                    }
                });
                this.app.reaction_page.invalidate(reaction);
            }
            if (this.create) {
                this.app.reaction_page.datatable.row.add(reaction);
                this.app.reaction_page.datatable.sort();
                this.app.reaction_page.datatable.draw();
                this.app.settings_page.refresh();
            }
            this.app.history_page.refresh();
            $(this.dialog_element)["modal"]("hide");
            return true;
        }
    }
    exports.Dialog = Dialog;
});
