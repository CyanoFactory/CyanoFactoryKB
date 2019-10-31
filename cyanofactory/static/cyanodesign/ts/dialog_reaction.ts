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
                        <div class="form-group reaction-substrates">
                            <label for="reaction-substrates">Substrates</label>
                        </div>
                        <div class="form-group reaction-products">
                            <label for="reaction-products">Products</label>
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

let template_metabolite = document.createElement('template');
template_metabolite.innerHTML = `
<div class="form-group form-inline metabolites">
<input class="form-control metabolite-amount" style="width:10%;" type="text">
<input class="form-control metabolite-id" style="width:85%" type="text">
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
    readonly substrates: HTMLDivElement;
    readonly products: HTMLDivElement;
    readonly enabled: ElementWrapper<boolean>;
    readonly reversible: ElementWrapper<boolean>;
    readonly constrained: ElementWrapper<boolean>;
    readonly constrained_min: ElementWrapper<string>;
    readonly constrained_max: ElementWrapper<string>;

    readonly obj_options: any = {
        maxItems: 1,
        valueField: 'id',
        searchField: ['name', 'id'],
        create: function (input: string) {
            let met = new mm.Metabolite();
            met.id = input; // todo: sluggify
            met.name = input;
            return met;
        },
        render: {
            item: function(item: mm.Reaction, escape: any) {
                return "<div>" + escape(item.get_name_or_id()) + "</div>";
            },
            option: function(item: mm.Reaction, escape: any) {
                return "<div>" + escape(item.get_name_or_id()) + "</div>";
            }
        },
        placeholder: "Select metabolite to add"
    };

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
        this.substrates = <HTMLDivElement>this.dialog_element.getElementsByClassName("reaction-substrates")[0]!;
        this.products = <HTMLDivElement>this.dialog_element.getElementsByClassName("reaction-products")[0]!;
        this.enabled = ElementWrapper.byClass<boolean>("enabled", this.dialog_element);
        this.reversible = ElementWrapper.byClass<boolean>("reversible", this.dialog_element);
        this.constrained = ElementWrapper.byClass<boolean>("constrained", this.dialog_element);
        this.constrained_min = ElementWrapper.byClass<string>("constrained-min", this.dialog_element);
        this.constrained_max = ElementWrapper.byClass<string>("constrained-max", this.dialog_element);
    }

    show(reaction: mm.Reaction) {
        if (reaction == null) {
            this.item = new mm.Reaction();
            reaction = this.item;
            this.create = true;
            this.item.makeUnconstrained(this.app.model);
        } else {
            this.item = reaction;
            this.create = false;
        }

        // Copy to dialog
        this.id.value = this.item.id;
        this.name.value = this.item.name;
        this.enabled.value = this.item.enabled;
        this.reversible.value = this.item.reversible;
        this.constrained.value = this.item.isConstrained(this.app.model);
        this.constrained_min.value = this.item.lower_bound.toString();
        this.constrained_max.value = this.item.upper_bound.toString();

        // cleanup
        $(this.dialog_element).find(".metabolites").remove();

        // Update dropdowns
        for (const metref of reaction.substrates) {
            this.addSubstrate(metref);
        }
        this.addSubstrate(null);

        for (const metref of reaction.products) {
            this.addProduct(metref);
        }
        this.addProduct(null);

        // constraint div visibility
        const constraint_div: JQuery<HTMLElement> = $(this.dialog_element).find(".constraints-min-max");

        if (reaction.isConstrained(this.app.model)) {
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

        for (const met of this.substrates.getElementsByClassName("metabolites")) {
            const amount = DialogHelper.getFloat(met.getElementsByClassName("metabolite-amount")[0]!.value);
            valid = valid && DialogHelper.checkBool(met.getElementsByClassName("metabolite-amount")[0]!, Number.isFinite(amount), "Value must be numerc");
        }
        for (const met of this.products.getElementsByClassName("metabolites")) {
            const amount = DialogHelper.getFloat(met.getElementsByClassName("metabolite-amount")[0]!.value);
            valid = valid && DialogHelper.checkBool(met.getElementsByClassName("metabolite-amount")[0]!, Number.isFinite(amount), "Value must be numerc");
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
            reaction.isConstrained(this.app.model) != constraints_enabled ||
            reaction.lower_bound != cmin_float ||
            reaction.upper_bound != cmax_float;

        let old_id = reaction.id;
        reaction.updateId(this.id.value, this.app.model);

        reaction.name = this.name.value;
        reaction.enabled = this.enabled.value;
        reaction.reversible = this.reversible.value;

        /* constraints */
        if (!constraints_enabled) {
            reaction.makeUnconstrained(this.app.model);
        } else {
            reaction.lower_bound = cmin_float;
            reaction.upper_bound = cmax_float;
        }

        /* substrates & products */

        // remove all metabolites
        reaction.substrates = [];
        reaction.products = [];

        // add metabolites
        const met_adder = (mets: any, target: mm.MetaboliteReference[]) => {
            for (const i of mets.getElementsByClassName("metabolites")) {
                const amount = i.getElementsByClassName("metabolite-amount")[0]!.value;
                const id = i.getElementsByClassName("metabolite-id")[0]!.value;

                if (id == "") {
                    continue;
                }

                let met: mm.Metabolite | null = this.app.model.metabolite.get("id", id);
                if (met == null) {
                    // is a new metabolite
                    met = new mm.Metabolite();
                    met.id = id; // Todo: Sluggify
                    met.name = id;
                    this.app.metabolite_page.datatable.row.add(<any>met);
                    this.app.model.metabolites.push(met);

                    this.app.history_manager.push({
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
                metref.id = id;
                metref.stoichiometry = amount;

                target.push(metref);
            }
        };

        met_adder(this.substrates, reaction.substrates);
        met_adder(this.products, reaction.products);

        this.app.metabolite_page.datatable.sort();
        this.app.metabolite_page.datatable.draw();

        reaction.updateMetaboliteReference(this.app.model);

        if (any_changed) {
            this.app.history_manager.push({
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

    private addSubstrate(metabolite: mm.MetaboliteReference | null) {
        let child: any = template_metabolite.content.cloneNode(true);
        const amount = child.children[0].children[0];
        const name = child.children[0].children[1];

        this.substrates.appendChild(child);

        amount.value = "1";

        $(name).selectize(this.obj_options);

        let selectize: any = name.selectize;
        for (const met of this.app.model.metabolites) {
            selectize.addOption(met);
        }

        if (metabolite != null) {
            amount.value = metabolite.stoichiometry;
            selectize.setValue(metabolite.id);
        }

        const self: Dialog = this;

        selectize.on("blur", function() {
            if (this.$dropdown.get()[0] == $(self.substrates).find(".metabolite-id").last().get()[0]) {
                if (this.getValue() != "") {
                    self.addSubstrate(null);
                }
            } else if (this.getValue() == "") {
                this.$dropdown.get()[0].parentNode.parentNode.remove();
            }
        });
    }

    private addProduct(metabolite: mm.MetaboliteReference | null) {
        let child: any = template_metabolite.content.cloneNode(true);
        const amount = child.children[0].children[0];
        const name = child.children[0].children[1];

        this.products.appendChild(child);

        amount.value = "1";

        $(name).selectize(this.obj_options);

        let selectize: any = name.selectize;
        for (const met of this.app.model.metabolites) {
            selectize.addOption(met);
        }

        if (metabolite != null) {
            amount.value = metabolite.stoichiometry;
            selectize.setValue(metabolite.id);
        }

        const self: Dialog = this;

        selectize.on("blur", function() {
            if (this.$dropdown.get()[0] == $(self.products).find(".metabolite-id").last().get()[0]) {
                if (this.getValue() != "") {
                    self.addProduct(null);
                }
            } else if (this.getValue() == "") {
                this.$dropdown.get()[0].parentNode.parentNode.remove();
            }
        });
    }
}
