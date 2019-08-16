import * as app from "./app"
import * as mm from "./metabolic_model";
import { DialogHelper, ElementWrapper } from "./dialog_helper";
import * as $ from "jquery";
import "datatables.net";

let template = document.createElement('template');
template.innerHTML = `
<div class="dialog-metabolite modal fade" tabindex="-1" role="dialog" aria-labelledby="dialog-add-metabolite-label"
     aria-hidden="true">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span
                        aria-hidden="true">&times;</span></button>
                <h4 class="modal-title" id="dialog-add-metabolite-label">Edit Metabolite</h4>
            </div>
            <div class="modal-body">
                <form>
                    <fieldset>
                        <div class="form-group">
                            <label for="metabolite-id">Identifier</label>
                            <input type="text" class="metabolite-id form-control">
                        </div>
                    
                        <div class="form-group">
                            <label for="add-metabolite-name">Name</label>
                            <input type="text" class="metabolite-name form-control">
                        </div>
                        
                        <div class="form-group">
                            <label class="control-label">Compartment</label>
                            <select class="metabolite-compartment" placeholder="Select a compartment"></select>                        
                        </div>

                        <div class="form-group">
                            <label for="add-chemical">Chemical formula</label>
                            <input type="text" class="metabolite-chemical form-control">
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
    item: mm.Metabolite = null;
    create: boolean = true;

    readonly id: ElementWrapper<string>;
    readonly name: ElementWrapper<string>;
    readonly compartment: ElementWrapper<string>;
    readonly formula: ElementWrapper<string>;

    constructor(app: app.AppManager) {
        this.app = app;

        document.body.appendChild(template.content.cloneNode(true));
        this.dialog_element = <HTMLElement>document.body.getElementsByClassName("dialog-metabolite")[0]!;
        const self: Dialog = this;
        $(this.dialog_element).find(".btn-primary").click(function () {
            self.validate();
        });

        this.id = ElementWrapper.byClass<string>("metabolite-id", this.dialog_element);
        this.name = ElementWrapper.byClass<string>("metabolite-name", this.dialog_element);
        this.compartment = ElementWrapper.byClass<string>("metabolite-compartment", this.dialog_element);
        this.formula = ElementWrapper.byClass<string>("metabolite-chemical", this.dialog_element);
    }

    show(metabolite: mm.Metabolite = null) {
        if (metabolite == null) {
            this.item = new mm.Metabolite();
            this.create = true;
        } else {
            this.item = metabolite;
            this.create = false;
        }

        // Copy to dialog
        this.id.value = this.item.id;
        this.name.value = this.item.name;
        this.formula.value = this.item.formula;

        // add compartments
        $(this.compartment.element).selectize({
            valueField: 'id',
            labelField: 'name',
            searchField: ['name', 'id'],
            options: this.app.model.compartments,
            persist: false,
            render: {
                item: function (item: mm.Compartment, escape: any) {
                    return "<div>" + escape(item.get_name_or_id()) + "</div>";
                },
                option: function (item: mm.Compartment, escape: any) {
                    return "<div>" + escape(item.get_name_or_id()) + "</div>";
                }
            }
        });
        let sel: any = this.compartment.element.selectize;
        sel.clearOptions();
        for (let c of this.app.model.compartments) {
            sel.addOption(c);
        }
        sel.clear();

        sel.addItem(this.item.compartment);

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
        let valid = DialogHelper.checkId(this.id.element);
        valid = valid && DialogHelper.checkLength(this.name.element, "Name", 1, 255);
        valid = valid && DialogHelper.checkRegexpPos(this.name.element, /(^[0-9])/, "Name must not begin with a number");
        valid = valid && DialogHelper.checkRegexpPos(this.name.element, /(^<?\-> | <?\-> | <?\->$)/, "Name must not contain lonely <-> or ->");
        valid = valid && DialogHelper.checkRegexpPos(this.name.element, /(^\+ | \+ )/, "Lonely + only allowed at end of name");

        let met_with_id = this.app.model.metabolite.get("id", this.id.value);

        if (met_with_id != null) {
            valid = valid && DialogHelper.checkBool(this.id.element, met_with_id == this.item, "Identifier already in use");
        }

        if (!valid) {
            return false;
        }

        let metabolite: mm.Metabolite = this.item;
        if (this.create) {
            metabolite.id = this.id.value;
            this.app.model.metabolites.push(metabolite);
        }

        let any_changed: boolean = this.create ||
            metabolite.id != this.id.value ||
            metabolite.name != this.name.value ||
            metabolite.compartment != this.compartment.value ||
            metabolite.formula != this.formula.value;

        let old_id = metabolite.id;
        metabolite.updateId(this.id.value, this.app.model);

        metabolite.compartment = this.compartment.value;
        metabolite.name = this.name.value;
        metabolite.formula = this.formula.value;

        if (any_changed) {
            this.app.command_list.push({
                "type": "metabolite",
                "op": this.create ? "add" : "edit",
                "id": old_id,
                "object": {
                    "id": metabolite.id,
                    "name": metabolite.name,
                    "compartment": metabolite.compartment,
                    "formula": metabolite.formula
                }
            });

            this.app.metabolite_page.invalidate(metabolite);
            this.app.compartment_page.invalidate();
        }

        if (this.create) {
            this.app.metabolite_page.datatable.row.add(metabolite);
            this.app.metabolite_page.datatable.sort();
            this.app.metabolite_page.datatable.draw();
        }

        this.app.history_page.refresh();

        $(this.dialog_element)["modal"]("hide");

        return true;
    }
}
