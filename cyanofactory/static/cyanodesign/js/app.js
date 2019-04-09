define(["require", "exports", "./page_reactions", "./page_metabolites", "./page_stoichiometry", "./page_settings", "./page_simulation", "./dialog_reaction", "./dialog_reaction_bulkadd", "./dialog_metabolite"], function (require, exports, reactions, metabolites, stoichiometry, settings, simulation, dialog_reaction, dialog_reaction_bulkadd, dialog_metabolite) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class AppManager {
        constructor(model) {
            this.command_list = [];
            this.dialog_reaction = new dialog_reaction.Dialog(model);
            this.dialog_reaction_bulk = new dialog_reaction_bulkadd.Dialog(model);
            this.dialog_metabolite = new dialog_metabolite.Dialog(model);
            this.reaction_page = new reactions.Page(document.getElementById("reaction-tab"), this);
            this.metabolite_page = new metabolites.Page(document.getElementById("metabolite-tab"), this);
            this.stoichiometry_page = new stoichiometry.Page(document.getElementById("chemical-tab"), this);
            this.settings_page = new settings.Page(document.getElementById("settings-tab"), this);
            this.simulation_page = new simulation.Page(document.getElementById("simulation-tab"), this);
            this.model = model;
        }
    }
    exports.AppManager = AppManager;
    function run(model) {
        app = new AppManager(model);
        app.reaction_page.update();
        app.metabolite_page.update();
        app.stoichiometry_page.update();
        app.settings_page.update();
        app.simulation_page.update();
    }
    exports.run = run;
});
