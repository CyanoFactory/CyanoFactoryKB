import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

export class Dialog {
    readonly app: app.AppManager;
    readonly model: mm.Model;
    readonly dialog_element: HTMLElement;
    private item: mm.Reaction = null;

    constructor(app: app.AppManager) {
        this.app = app;

        this.dialog_element = <HTMLElement>document.getElementById("dialog-saveas-model");
        const self: Dialog = this;
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

    ok_clicked(): void {
        this.app.request_handler.saveas(
            <string>$(this.dialog_element).find("#id_saveas_summary").val(),
            <string>$(this.dialog_element).find("#id_saveas_name").val());
        $(this.dialog_element)["modal"]("hide");
    }
}
