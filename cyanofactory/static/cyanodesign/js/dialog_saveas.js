define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class Dialog {
        constructor(app) {
            this.item = null;
            this.app = app;
            this.dialog_element = document.getElementById("dialog-saveas-model");
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
            this.app.request_handler.saveas($(this.dialog_element).find("#id_saveas_summary").val(), $(this.dialog_element).find("#id_saveas_name").val());
            $(this.dialog_element)["modal"]("hide");
        }
    }
    exports.Dialog = Dialog;
});
