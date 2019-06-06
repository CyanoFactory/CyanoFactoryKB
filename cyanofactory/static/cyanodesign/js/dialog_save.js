define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class Dialog {
        constructor(app) {
            this.item = null;
            this.app = app;
            this.dialog_element = document.getElementById("dialog-save-model");
            const self = this;
            $(this.dialog_element).find(".btn-primary").click(function () {
                self.ok_clicked();
            });
        }
        show() {
            // cleanup
            $(this.dialog_element).find(".modal-content");
            // FIXME
            $("#id_save_summary").val("");
            // show
            $(this.dialog_element)["modal"]("show");
        }
        ok_clicked() {
            const form = $("#dialog-save-model").find("form");
            this.app.request_handler.save(form.find("#id_save_summary").val());
            $(this.dialog_element)["modal"]("hide");
        }
    }
    exports.Dialog = Dialog;
});
