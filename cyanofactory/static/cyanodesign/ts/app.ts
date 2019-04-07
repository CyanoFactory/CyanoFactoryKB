import * as mm from "./metabolic_model"
import * as $ from "jquery";
import * as reactions from "./page_reactions"
import * as metabolites from "./page_metabolites"
import * as stoichiometry from "./page_stoichiometry"
import * as settings from "./page_settings"
import * as simulation from "./page_simulation"
import * as dialog_reaction from "./dialog_reaction"
import * as dialog_reaction_bulkadd from "./dialog_reaction_bulkadd"
import * as dialog_metabolite from "./dialog_metabolite"

export class AppManager {
    readonly dialog_reaction: dialog_reaction.Dialog;
    readonly dialog_reaction_bulk: dialog_reaction_bulkadd.Dialog;
    readonly dialog_metabolite: dialog_metabolite.Dialog;
    readonly model: mm.Model;

    constructor(model: mm.Model) {
        this.dialog_reaction = new dialog_reaction.Dialog(model);
        this.dialog_reaction_bulk = new dialog_reaction_bulkadd.Dialog(model);
        this.dialog_metabolite = new dialog_metabolite.Dialog(model);
        this.model = model;
    }
}

export function run(model: mm.Model) {
    let app = new AppManager(model);

    let reaction_page = new reactions.Page(document.getElementById("reaction-tab")!, app);
    reaction_page.update();

    let metabolite_page = new metabolites.Page(document.getElementById("metabolite-tab")!, app);
    metabolite_page.update();

    let stoichiometry_page = new stoichiometry.Page(document.getElementById("chemical-tab")!, app);
    stoichiometry_page.update();

    let settings_page = new settings.Page(document.getElementById("settings-tab")!, app);
    settings_page.update();

    let simulation_page = new simulation.Page(document.getElementById("simulation-tab")!, app);
    simulation_page.update();
}
