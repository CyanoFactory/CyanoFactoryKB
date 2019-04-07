define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
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
                <h4 class="modal-title" id="dialog-add-metabolite-label">Create Metabolite</h4>
            </div>
            <div class="modal-body">
                <form>
                    <fieldset>
                        <div class="form-group">
                            <label for="add-metabolite-name">Name</label>
                            <input type="text" id="add-metabolite-name" class="form-control">
                        </div>

                        <div class="checkbox">
                            <input type="checkbox" id="add-external">
                            <label for="add-external">External</label>
                        </div>

                        <div class="form-group">
                            <label for="add-chemical">Chemical formula</label>
                            <input type="text" id="metabolite-chemical" class="form-control">
                        </div>
                    </fieldset>
                </form>
            </div>
            <div class="modal-footer">
                <button type="button" class="btn btn-primary">Create</button>
                <button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
            </div>
        </div>
    </div>
</div>
`;
    class Dialog {
        constructor(model) {
            document.body.appendChild(template.content.cloneNode(true));
            this.dialog_element = document.body.getElementsByClassName("dialog-metabolite")[0];
            this.model = model;
        }
        show() {
            /*
            function showAddMetaboliteDialog() {
                dialog_metabolite_add.find("div").removeClass("has-error");
                dialog_metabolite_add.find(".help-block").remove();
                dialog_metabolite_add.find("#add-metabolite-name").val("New Metabolite");
                dialog_metabolite_add.find("#add-external").prop("checked", false);
                dialog_metabolite_add.modal("show");
            }
    
            function showEditMetaboliteDialog(item) {
                dialog_metabolite_edit.find("div").removeClass("has-error");
                dialog_metabolite_edit.find(".help-block").remove();
                dialog_metabolite_edit.find("#metabolite-name").val(item.name);
                dialog_metabolite_edit.find("#external").prop("checked", model.metabolite.get("id", item.id).isExternal());
                dialog_metabolite_edit.find("#metabolite-chemical").val(item.chemical);
                dialog_metabolite_edit.data("object", item);
                dialog_metabolite_edit.modal("show");
            }*/
            $(this.dialog_element)["modal"]("show");
        }
    }
    exports.Dialog = Dialog;
});
