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
                            <label for="add-metabolite-name">Identifier</label>
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
        }
        get id() {
            return $(this.dialog_element).find(".metabolite-id").val();
        }
        set id(id) {
            $(this.dialog_element).find(".metabolite-id").val(id);
        }
        get name() {
            return $(this.dialog_element).find(".metabolite-name").val();
        }
        set name(name) {
            $(this.dialog_element).find(".metabolite-name").val(name);
        }
        get external() {
            return $(this.dialog_element).find(".metabolite-external").prop("checked");
        }
        set external(external) {
            $(this.dialog_element).find(".metabolite-external").prop("checked", external);
        }
        get formula() {
            return $(this.dialog_element).find(".metabolite-chemical").val();
        }
        set formula(name) {
            $(this.dialog_element).find(".metabolite-chemical").val(name);
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
            this.id = this.item.id;
            this.name = this.item.name;
            this.external = this.item.isExternal();
            this.formula = this.item.formula;
            // clean up
            $(this.dialog_element).find("div").removeClass("has-error");
            $(this.dialog_element).find(".help-block").remove();
            // show
            $(this.dialog_element)["modal"]("show");
        }
        validate() {
            let name = $(this.dialog_element).find(".metabolite-name");
            let valid = dialog_helper_1.DialogHelper.checkLength(name, "Name", 1, 100);
            valid = valid || dialog_helper_1.DialogHelper.checkRegexpPos(name, /(^[0-9])/, "Name must not begin with a number");
            valid = valid || dialog_helper_1.DialogHelper.checkRegexpPos(name, /(^<?\-> | <?\-> | <?\->$)/, "Name must not contain lonely <-> or ->");
            valid = valid || dialog_helper_1.DialogHelper.checkRegexpPos(name, /(^\+ | \+ )/, "Lonely + only allowed at end of name");
            let met_with_id = this.model.metabolite.get("id", this.id);
            if (met_with_id != null) {
                valid = valid || dialog_helper_1.DialogHelper.checkBool(name, met_with_id == this.item, "Name already in use");
            }
            if (valid) {
                var metabolite = dialog_metabolite_edit.data("object");
                var old_name = metabolite.id;
                metabolite.updateName(name.val());
                var old_external = metabolite.external;
                metabolite.external = external.prop("checked");
                if (old_name != metabolite.name || old_external != metabolite.external) {
                    command_list.push({
                        "type": "metabolite",
                        "op": "edit",
                        "id": old_name,
                        "object": {
                            "id": metabolite.id,
                            "name": metabolite.name,
                            "external": metabolite.external,
                            "chemical": metabolite.chemical
                        }
                    });
                }
                metabolite.invalidate();
                dialog_metabolite_edit.modal("hide");
            }
            /*
            
                    var name = dialog_metabolite_add.find("#add-metabolite-name");
                    var external = dialog_metabolite_add.find("#add-external");
            
                    var valid = checkLength(name, "Name", 1, 100);
                    valid &= checkRegexpPos(name, /(^[0-9])/, "Name must not begin with a number");
                    valid &= checkRegexpPos(name, /(^<?\-> | <?\-> | <?\->$)/, "Name must not contain lonely <-> or ->");
                    valid &= checkRegexpPos(name, /(^\+ | \+ )/, "Lonely + only allowed at end of name");
            
                    var met_idx = Metabolite.indexByName(name.val());
            
                    if (met_idx != -1) {
                        valid = checkBool(name, false, "Name already in use");
                    }
            
                    if (valid) {
                        var metabolite = new Metabolite(name.val(), name.val(), external.prop("checked"));
                        model.metabolites.push(metabolite);
            
                        command_list.push({
                            "type": "metabolite",
                            "op": "add",
                            "id": metabolite.id,
                            "object": {
                                "name": metabolite.name,
                                "external": metabolite.external,
                                "chemical": metabolite.chemical
                            }
                        });
            
                        datatable_metabolites.row.add(metabolite);
                        datatable_metabolites.sort();
                        datatable_metabolites.draw();
            
                        dialog_metabolite_add.modal("hide");
                    }*/
        }
    }
    exports.Dialog = Dialog;
});
