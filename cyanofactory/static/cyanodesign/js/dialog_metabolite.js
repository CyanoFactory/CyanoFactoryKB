define(["require", "exports", "./metabolic_model", "./dialog_helper", "jquery", "datatables.net"], function (require, exports, mm, dialog_helper_1, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
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

                        <div class="checkbox">
                            <input type="checkbox" class="metabolite-external">
                            <label for="add-external">External</label>
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
    class Dialog {
        constructor(model) {
            this.item = null;
            this.create = true;
            document.body.appendChild(template.content.cloneNode(true));
            this.dialog_element = document.body.getElementsByClassName("dialog-metabolite")[0];
            this.model = model;
            const self = this;
            $(this.dialog_element).find(".btn-primary").click(function () {
                self.validate();
            });
            this.id = dialog_helper_1.ElementWrapper.byClass("metabolite-id", this.dialog_element);
            this.name = dialog_helper_1.ElementWrapper.byClass("metabolite-name", this.dialog_element);
            this.external = dialog_helper_1.ElementWrapper.byClass("metabolite-external", this.dialog_element);
            this.formula = dialog_helper_1.ElementWrapper.byClass("metabolite-chemical", this.dialog_element);
        }
        show(metabolite = null) {
            if (metabolite == null) {
                this.item = new mm.Metabolite();
                this.create = true;
            }
            else {
                this.item = metabolite;
                this.create = false;
            }
            // Copy to dialog
            this.id.value = this.item.id;
            this.name.value = this.item.name;
            this.external.value = this.item.isExternal();
            this.formula.value = this.item.formula;
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
            let valid = dialog_helper_1.DialogHelper.checkId(this.id.element);
            valid = valid && dialog_helper_1.DialogHelper.checkLength(this.name.element, "Name", 1, 255);
            valid = valid && dialog_helper_1.DialogHelper.checkRegexpPos(this.name.element, /(^[0-9])/, "Name must not begin with a number");
            valid = valid && dialog_helper_1.DialogHelper.checkRegexpPos(this.name.element, /(^<?\-> | <?\-> | <?\->$)/, "Name must not contain lonely <-> or ->");
            valid = valid && dialog_helper_1.DialogHelper.checkRegexpPos(this.name.element, /(^\+ | \+ )/, "Lonely + only allowed at end of name");
            let met_with_id = this.model.metabolite.get("id", this.id.value);
            if (met_with_id != null) {
                valid = valid && dialog_helper_1.DialogHelper.checkBool(this.id.element, met_with_id == this.item, "Identifier already in use");
            }
            if (!valid) {
                return false;
            }
            let metabolite = this.item;
            let any_changed = metabolite.id != this.id.value ||
                metabolite.name != this.name.value ||
                metabolite.isExternal() != this.external.value ||
                metabolite.formula != this.formula.value;
            let old_id = metabolite.id;
            metabolite.updateId(this.id.value, this.model);
            // FIXME: Compartment configuration?
            metabolite.compartment = this.external.value ? "e" : "c";
            metabolite.name = this.name.value;
            metabolite.formula = this.formula.value;
            if (any_changed) {
                app.command_list.push({
                    "type": "metabolite",
                    "op": "edit",
                    "id": old_id,
                    "object": {
                        "id": metabolite.id,
                        "name": metabolite.name,
                        "external": metabolite.external,
                        "formula": metabolite.formula
                    }
                });
                app.metabolite_page.invalidate(metabolite);
            }
            $(this.dialog_element)["modal"]("hide");
            return true;
        }
    }
    exports.Dialog = Dialog;
});
