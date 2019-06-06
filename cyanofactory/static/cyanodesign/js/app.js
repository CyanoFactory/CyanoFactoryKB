define(["require", "exports", "jquery", "./page_reactions", "./page_metabolites", "./page_settings", "./page_simulation", "./page_history", "./dialog_reaction", "./dialog_reaction_bulkadd", "./dialog_reaction_delete", "./dialog_metabolite", "./dialog_save", "./request_handler"], function (require, exports, $, reactions, metabolites, settings, simulation, history, dialog_reaction, dialog_reaction_bulkadd, dialog_reaction_delete, dialog_metabolite, dialog_save, request_handler_1) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class AppManager {
        constructor(mm_cls, model_json, urls) {
            this.command_list = [];
            this.model = new mm_cls.Model();
            this.model.fromJson(model_json);
            this.original_model = new mm_cls.Model();
            this.original_model.fromJson(model_json);
            this.urls = urls;
            this.request_handler = new request_handler_1.RequestHandler(this, -1);
            this.dialog_reaction = new dialog_reaction.Dialog(this);
            this.dialog_reaction_bulk = new dialog_reaction_bulkadd.Dialog(this);
            this.dialog_reaction_delete = new dialog_reaction_delete.Dialog(this);
            this.dialog_metabolite = new dialog_metabolite.Dialog(this);
            this.dialog_save = new dialog_save.Dialog(this);
            this.reaction_page = new reactions.Page(document.getElementById("reaction-tab"), this);
            this.metabolite_page = new metabolites.Page(document.getElementById("metabolite-tab"), this);
            //this.stoichiometry_page = new stoichiometry.Page(document.getElementById("chemical-tab")!, this);
            this.settings_page = new settings.Page(document.getElementById("settings-tab"), this);
            this.history_page = new history.Page(document.getElementById("history-tab"), this);
            this.simulation_page = new simulation.Page(document.getElementById("simulation-tab"), this);
        }
    }
    exports.AppManager = AppManager;
    function run(mm_cls, urls) {
        $.ajax({
            url: urls.get_reactions,
            context: document.body
        }).done(function (x) {
            app = new AppManager(mm_cls, x, urls);
            app.reaction_page.init();
            app.metabolite_page.init();
            //app.stoichiometry_page.init();
            app.settings_page.init();
            app.history_page.init();
            app.simulation_page.init();
            $(".create-enzyme-button").on("click", function () {
                app.dialog_reaction.show(null);
            });
            $(".create-metabolite-button").on("click", function () {
                app.dialog_metabolite.show(null);
            });
            $(".create-enzyme-bulk-button").on("click", function () {
                app.dialog_reaction_bulk.show();
            });
            $("#design-save").on("click", function () {
                app.dialog_save.show();
            });
            app.request_handler.endRequest($("#content"));
        });
    }
    exports.run = run;
});
