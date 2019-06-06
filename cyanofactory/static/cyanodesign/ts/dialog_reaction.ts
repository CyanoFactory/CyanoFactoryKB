import * as app from "./app"
import * as mm from "./metabolic_model";
import { DialogHelper, ElementWrapper } from "./dialog_helper";
import * as $ from "jquery";
import "datatables.net";
import "selectize";

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
    readonly app: app.AppManager;
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

    public substrate_refs: mm.MetaboliteReference[] = [];
    public product_refs: mm.MetaboliteReference[] = [];

    constructor(app: app.AppManager) {
        this.app = app;

        document.body.appendChild(template.content.cloneNode(true));
        this.dialog_element = <HTMLElement>document.body.getElementsByClassName("dialog-reaction")[0]!;
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

    show(reaction: mm.Reaction) {
        if (reaction == null) {
            this.item = new mm.Reaction();
            reaction = this.item;
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

        this.substrate_refs = [];
        this.product_refs = [];

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

        // FIXME: HACK
        for (let substrate of this.item.substrates) {
            for (let ref of this.substrate_refs) {
                if (substrate.id == ref.id) {
                    ref.stoichiometry = substrate.stoichiometry;
                    break;
                }
            }
        }

        for (let product of this.item.products) {
            for (let ref of this.product_refs) {
                if (product.id == ref.id) {
                    ref.stoichiometry = product.stoichiometry;
                    break;
                }
            }
        }

        // Update dropdowns
        let obj_options: any = {
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
                },
                'add_option_hook' : {}
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

        $(this.substrates.element).selectize(obj_options);
        let substrates_selectize: any = this.substrates.element.selectize;
        substrates_selectize.clearOptions();
        for (let ref of this.substrate_refs) {
            substrates_selectize.addOption(ref);
        }

        substrates_selectize.clear();

        for (const metref of reaction.substrates) {
            substrates_selectize.addItem(metref.id);
        }

        obj_options["options"] = this.product_refs;
        $(this.products.element)["selectize"](obj_options);
        let products_selectize: any = this.products.element.selectize;
        products_selectize.clearOptions();
        for (let ref of this.product_refs) {
            products_selectize.addOption(ref);
        }

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

        let reac_with_id = this.app.model.reaction.get("id", this.id.value);

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
            this.app.model.reactions.push(reaction);
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
        reaction.updateId(this.id.value, this.app.model);

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
        let substrates_selectize: any = this.substrates.element.selectize;
        let products_selectize: any = this.products.element.selectize;

        // remove all metabolites
        reaction.substrates = [];
        reaction.products = [];

        // add metabolites
        const met_adder = (mets: any, target: mm.MetaboliteReference[]) => {
            for (const i of mets.items) {
                let met: mm.Metabolite | null = this.app.model.metabolite.get("id", i);
                if (met == null) {
                    // is a new metabolite
                    met = new mm.Metabolite();
                    met.id = i; // Todo: Sluggify
                    met.name = i;
                    this.app.metabolite_page.datatable.row.add(<any>met);
                    this.app.model.metabolites.push(met);

                    this.app.command_list.push({
                        "type": "metabolite",
                        "op": "add",
                        "id": met.id,
                        "object": {
                            "id": met.id,
                            "name": met.name,
                            "external": met.external,
                            "formula": met.formula
                        }
                    });
                }

                let metref = new mm.MetaboliteReference();
                metref.id = i;
                metref.stoichiometry = mets.options[metref.id].stoichiometry;

                target.push(metref);
            }
        };

        met_adder(substrates_selectize, reaction.substrates);
        met_adder(products_selectize, reaction.products);

        this.app.metabolite_page.datatable.sort();
        this.app.metabolite_page.datatable.draw();

        reaction.updateMetaboliteReference(this.app.model);

        if (any_changed) {
            this.app.command_list.push({
                "type": "reaction",
                "op": this.create ? "add" : "edit",
                "id": old_id,
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
