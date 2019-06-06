define(["require", "exports", "jquery"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class RequestHandler {
        constructor(app, revision) {
            this.app = app;
            this.revision = revision;
        }
        beginRequest(show_on = null) {
            if (show_on == null)
                waitIndicator.show($("#content"));
            else
                waitIndicator.show(show_on);
            let result = {
                "changes": JSON.stringify(this.app.command_list),
                "objectives": JSON.stringify([{
                        "id": this.app.settings_page.getObjective(),
                        "maximize": this.app.settings_page.maximizeObjective()
                    }]),
                "design_objectives": JSON.stringify([{
                        "id": this.app.settings_page.getDesignObjective(),
                        "maximize": true
                    }]),
                "target_reactions": JSON.stringify([{
                        "id": this.app.settings_page.getTargetObjective(),
                        "maximize": true
                    }]),
                "type": JSON.stringify(this.app.settings_page.getSimulationType()),
                //"display": JSON.stringify(design_objective_visible_combobox[0].selectize.getValue().split("\x00")),
                //"auto_flux": JSON.stringify($("#auto_flux").prop("checked"))
                "display": JSON.stringify([]),
                "auto_flux": JSON.stringify(true)
            };
            if (this.revision != -1) {
                result["revision"] = this.revision;
            }
            return result;
        }
        endRequest(hide_on = null) {
            if (hide_on == null)
                waitIndicator.hide($("#content"));
            else
                waitIndicator.hide(hide_on);
        }
        simulate() {
            let data = this.beginRequest();
            const self = this;
            $.ajax({
                url: this.app.urls.simulate,
                type: "POST",
                data: data
            }).done(function (x) {
                self.app.simulation_page.simulate(x);
                self.endRequest();
            }).fail(function (x) {
                //$("#visual_graph").hide();
                //$("#visual_fba").hide();
                self.app.simulation_page.notifyError(x.responseText);
                self.endRequest();
            });
        }
        save(summary) {
            let data = this.beginRequest();
            const self = this;
            data["summary"] = summary;
            $.ajax({
                url: self.app.urls.save,
                type: "POST",
                data: data
            }).done(function (x) {
                /*if (revision !== undefined) {
                    location.replace(self.app.urls.design);
                }*/
                /*$("#dialog-save-model")["modal"]("hide");
                endRequest($("#dialog-save-model").find(".modal-content"));
                $("#simulation-result").html("")*/
                self.app.command_list = [];
                self.endRequest();
            }).fail(function (x) {
                /*$("#dialog-save-model")["modal"]("hide");
                $("#visual_graph").hide();
                $("#visual_fba").hide();*/
                //notifyError(JSON.parse(x.responseText)["message"]);
                self.endRequest();
            });
        }
        saveas() {
            let form = $("#dialog-saveas-model").find("form");
            let data = this.beginRequest($("#dialog-saveas-model").find(".modal-content"));
            data["saveas_name"] = form.find("#id_saveas_name").val();
            data["saveas_summary"] = form.find("#id_saveas_summary").val();
            const self = this;
            $.ajax({
                url: this.app.urls.saveas,
                type: "POST",
                data: data
            }).done(function (x) {
                if (!(x['success'])) {
                    form.replaceWith(x['form_html']);
                }
                else {
                    location.replace(x['url']);
                }
                self.endRequest($("#dialog-saveas-model").find(".modal-content"));
            }).fail(function (x) {
                $("#visual_graph").hide();
                $("#visual_fba").hide();
                //DialogHelper.notnotifyError(x.responseText);
                self.endRequest($("#dialog-saveas-model").find(".modal-content"));
            });
        }
    }
    exports.RequestHandler = RequestHandler;
});
