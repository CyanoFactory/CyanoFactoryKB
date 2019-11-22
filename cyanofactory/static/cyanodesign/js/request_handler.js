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
                "changes": JSON.stringify(this.app.history_manager.current())
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
            }).done((x) => {
                $.ajax({
                    url: this.app.urls.get_revisions,
                }).done((x) => this.app.history_page.init(x));
                self.endRequest();
            }).fail(function (x) {
                self.endRequest();
            });
        }
        saveas(summary, new_id) {
            let form = $("#dialog-saveas-model").find("form");
            let data = this.beginRequest();
            data["saveas_name"] = new_id;
            data["saveas_summary"] = summary;
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
                self.endRequest($("#dialog-saveas-model").find(".modal-content"));
            });
        }
    }
    exports.RequestHandler = RequestHandler;
});
