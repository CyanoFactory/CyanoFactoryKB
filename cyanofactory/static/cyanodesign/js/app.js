define(["require", "exports", "jquery", "./page_reactions", "./page_metabolites", "./page_stoichiometry", "./page_settings", "./page_simulation", "./dialog_reaction", "./dialog_reaction_bulkadd", "./dialog_reaction_delete", "./dialog_metabolite", "./request_handler"], function (require, exports, $, reactions, metabolites, stoichiometry, settings, simulation, dialog_reaction, dialog_reaction_bulkadd, dialog_reaction_delete, dialog_metabolite, request_handler_1) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class AppManager {
        constructor(model, urls) {
            this.command_list = [];
            this.model = model;
            this.urls = urls;
            this.request_handler = new request_handler_1.RequestHandler(this, -1);
            this.dialog_reaction = new dialog_reaction.Dialog(this);
            this.dialog_reaction_bulk = new dialog_reaction_bulkadd.Dialog(this);
            this.dialog_reaction_delete = new dialog_reaction_delete.Dialog(this);
            this.dialog_metabolite = new dialog_metabolite.Dialog(this);
            this.reaction_page = new reactions.Page(document.getElementById("reaction-tab"), this);
            this.metabolite_page = new metabolites.Page(document.getElementById("metabolite-tab"), this);
            this.stoichiometry_page = new stoichiometry.Page(document.getElementById("chemical-tab"), this);
            this.settings_page = new settings.Page(document.getElementById("settings-tab"), this);
            this.simulation_page = new simulation.Page(document.getElementById("simulation-tab"), this);
        }
    }
    exports.AppManager = AppManager;
    function run(mm_cls, urls) {
        $.ajax({
            url: urls.get_reactions,
            context: document.body
        }).done(function (x) {
            let model = new mm_cls.Model();
            model.fromJson(x);
            app = new AppManager(model, urls);
            app.reaction_page.init();
            app.metabolite_page.init();
            app.stoichiometry_page.init();
            app.settings_page.init();
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
            app.request_handler.endRequest($("#content"));
        });
    }
    exports.run = run;
});
