import * as app from "./app"
import * as mm from "./metabolic_model";
import { DialogHelper, ElementWrapper } from "./dialog_helper";
import * as $ from "jquery";
import "datatables.net";

let template = document.createElement('template');
template.innerHTML = `
<div class="dialog-compartment modal fade" tabindex="-1" role="dialog" aria-labelledby="dialog-add-metabolite-label"
     aria-hidden="true">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span
                        aria-hidden="true">&times;</span></button>
                <h4 class="modal-title" id="dialog-add-compartment-label">Edit Compartment</h4>
            </div>
            <div class="modal-body">
                <form>
                    <fieldset>
                        <div class="form-group">
                            <label for="compartment-id">Identifier</label>
                            <input type="text" class="compartment-id form-control">
                        </div>
                    
                        <div class="form-group">
                            <label for="compartment-name">Name</label>
                            <input type="text" class="compartment-name form-control">
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
    item: mm.Compartment = null;
    create: boolean = true;

    readonly id: ElementWrapper<string>;
    readonly name: ElementWrapper<string>;

    constructor(app: app.AppManager) {
        this.app = app;

        document.body.appendChild(template.content.cloneNode(true));
        this.dialog_element = <HTMLElement>document.body.getElementsByClassName("dialog-compartment")[0]!;
        const self: Dialog = this;
        $(this.dialog_element).find(".btn-primary").click(function () {
            self.validate();
        });

        this.id = ElementWrapper.byClass<string>("compartment-id", this.dialog_element);
        this.name = ElementWrapper.byClass<string>("compartment-name", this.dialog_element);
    }

    show(compartment: mm.Compartment = null) {
        if (compartment == null) {
            this.item = new mm.Compartment();
            this.create = true;
        } else {
            this.item = compartment;
            this.create = false;
        }

        // Copy to dialog
        this.id.value = this.item.id;
        this.name.value = this.item.name;

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

        let c_with_id = this.app.model.compartment.get("id", this.id.value);

        if (c_with_id != null) {
            valid = valid && DialogHelper.checkBool(this.id.element, c_with_id == this.item, "Identifier already in use");
        }

        if (!valid) {
            return false;
        }

        let compartment: mm.Compartment = this.item;
        if (this.create) {
            compartment.id = this.id.value;
            this.app.model.compartments.push(compartment);
        }

        let any_changed: boolean = this.create ||
            compartment.id != this.id.value ||
            compartment.name != this.name.value;

        let old_id = compartment.id;
        compartment.updateId(this.id.value, this.app.model);

        compartment.name = this.name.value;

        if (any_changed) {
            this.app.command_list.push({
                "type": "compartment",
                "op": this.create ? "add" : "edit",
                "id": old_id,
                "object": {
                    "id": compartment.id,
                    "name": compartment.name
                }
            });

            this.app.compartment_page.invalidate();
        }

        if (this.create) {
            this.app.compartment_page.datatable.row.add(compartment);
            this.app.compartment_page.datatable.sort();
            this.app.compartment_page.datatable.draw();
        }

        this.app.history_page.refresh();

        $(this.dialog_element)["modal"]("hide");

        return true;
    }
}
