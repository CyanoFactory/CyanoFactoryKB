define(["require", "exports", "./page_reactions", "./page_metabolites", "./page_stoichiometry", "./page_settings", "./page_simulation", "./dialog_reaction", "./dialog_reaction_bulkadd", "./dialog_metabolite"], function (require, exports, reactions, metabolites, stoichiometry, settings, simulation, dialog_reaction, dialog_reaction_bulkadd, dialog_metabolite) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class AppManager {
        constructor(model) {
            this.dialog_reaction = new dialog_reaction.Dialog(model);
            this.dialog_reaction_bulk = new dialog_reaction_bulkadd.Dialog(model);
            this.dialog_metabolite = new dialog_metabolite.Dialog(model);
            this.model = model;
        }
    }
    exports.AppManager = AppManager;
    function run(model) {
        let app = new AppManager(model);
        let reaction_page = new reactions.Page(document.getElementById("reaction-tab"), app);
        reaction_page.update();
        let metabolite_page = new metabolites.Page(document.getElementById("metabolite-tab"), app);
        metabolite_page.update();
        let stoichiometry_page = new stoichiometry.Page(document.getElementById("chemical-tab"), app);
        stoichiometry_page.update();
        let settings_page = new settings.Page(document.getElementById("settings-tab"), app);
        settings_page.update();
        let simulation_page = new simulation.Page(document.getElementById("simulation-tab"), app);
        simulation_page.update();
    }
    exports.run = run;
});
