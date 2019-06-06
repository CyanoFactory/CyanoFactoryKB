define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    let template = document.createElement('template');
    template.innerHTML = `
<div class="dialog-delete modal fade" tabindex="-1" role="dialog" aria-labelledby="dialog-delete-label"
     aria-hidden="true">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span
                        aria-hidden="true">&times;</span></button>
                <h4 class="modal-title dialog-delete-label">Delete <span class="reaction-name"></span></h4>
            </div>
            <div class="modal-body">
                Do you really want to delete reaction <span class="reaction-name"></span>?

                <div class="form-group">
                <div class="checkbox">
                    <input class="delete-unused-metabolites" type="checkbox">
                    <label for="dialog-delete-unused-metabolites">Delete unused metabolites</label>
                </div>
                </div>
            </div>
            <div class="modal-footer">
                <button type="button" class="btn btn-primary btn-danger">Delete</button>
                <button type="button" class="btn btn-default" data-dismiss="modal">Cancel</button>
            </div>
        </div>
    </div>
</div>
`;
    class Dialog {
        constructor(app) {
            this.item = null;
            this.app = app;
            document.body.appendChild(template.content.cloneNode(true));
            this.dialog_element = document.body.getElementsByClassName("dialog-delete")[0];
            const self = this;
            $(this.dialog_element).find(".btn-primary").click(function () {
                self.ok_clicked();
            });
        }
        show(reaction) {
            // cleanup
            $(this.dialog_element).find(".delete-unused-metabolites").prop("checked", false);
            this.item = reaction;
            $(this.dialog_element).find(".reaction-name").text(reaction.get_name_or_id());
            // show
            $(this.dialog_element)["modal"]("show");
        }
        ok_clicked() {
            let rmets = this.item.getMetabolites(this.app.model);
            let del_mets = $(this.dialog_element).find(".delete-unused-metabolites").prop("checked");
            this.app.reaction_page.datatable.row(this.app.model.reaction.checked_index("id", this.item.id)).remove();
            this.item.remove(this.app.model);
            for (let met of this.item.getMetabolites(this.app.model)) {
                this.app.metabolite_page.invalidate(met);
            }
            this.app.command_list.push({
                "op": "delete",
                "type": "reaction",
                "id": this.item.id,
                "object": {
                    "id": this.item.id
                }
            });
            if (del_mets) {
                rmets.filter((m) => {
                    // this.app.metabolite_page.invalidate(m);
                    return del_mets && m.isUnused();
                }).forEach((m) => {
                    this.app.metabolite_page.datatable.row(this.app.model.metabolite.checked_index("id", m.id)).remove();
                    m.remove(this.app.model);
                    this.app.command_list.push({
                        "type": "metabolite",
                        "op": "delete",
                        "id": m.id,
                        "object": {}
                    });
                });
            }
            this.app.reaction_page.refresh();
            this.app.metabolite_page.refresh();
            this.app.settings_page.refresh([this.item]);
            this.app.history_page.refresh();
            $(this.dialog_element)["modal"]("hide");
        }
    }
    exports.Dialog = Dialog;
});
