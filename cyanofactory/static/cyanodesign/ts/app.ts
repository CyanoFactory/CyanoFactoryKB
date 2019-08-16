import * as mm from "./metabolic_model"
import * as $ from "jquery";
import * as reactions from "./page_reactions"
import * as metabolites from "./page_metabolites"
import * as compartments from "./page_compartments"
import * as stoichiometry from "./page_stoichiometry"
import * as settings from "./page_settings"
import * as simulation from "./page_simulation"
import * as history from "./page_history"

import * as dialog_reaction from "./dialog_reaction"
import * as dialog_reaction_bulkadd from "./dialog_reaction_bulkadd"
import * as dialog_reaction_delete from "./dialog_reaction_delete"
import * as dialog_metabolite from "./dialog_metabolite"
import * as dialog_save from "./dialog_save"
import * as dialog_compartment from "./dialog_compartment"

import { RequestHandler } from "./request_handler"

declare var app : AppManager;

export class AppManager {
    readonly dialog_reaction: dialog_reaction.Dialog;
    readonly dialog_reaction_bulk: dialog_reaction_bulkadd.Dialog;
    readonly dialog_reaction_delete: dialog_reaction_delete.Dialog;
    readonly dialog_metabolite: dialog_metabolite.Dialog;
    readonly dialog_save: dialog_save.Dialog;
    readonly dialog_compartment: dialog_compartment.Dialog;
    readonly reaction_page: reactions.Page;
    readonly metabolite_page: metabolites.Page;
    readonly compartment_page: compartments.Page;
    readonly stoichiometry_page: stoichiometry.Page;
    readonly settings_page: settings.Page;
    readonly simulation_page: simulation.Page;
    readonly history_page: history.Page;
    readonly model: mm.Model;
    readonly original_model: mm.Model;
    readonly urls: any;
    readonly request_handler: RequestHandler;
    public command_list: any = [];

    constructor(mm_cls: any, model_json: any, urls: any) {
        this.model = new mm_cls.Model();
        this.model.fromJson(model_json);
        this.original_model = new mm_cls.Model();
        this.original_model.fromJson(model_json);

        this.urls = urls;
        this.request_handler = new RequestHandler(this, -1);

        this.dialog_reaction = new dialog_reaction.Dialog(this);
        this.dialog_reaction_bulk = new dialog_reaction_bulkadd.Dialog(this);
        this.dialog_reaction_delete = new dialog_reaction_delete.Dialog(this);
        this.dialog_metabolite = new dialog_metabolite.Dialog(this);
        this.dialog_save = new dialog_save.Dialog(this);
        this.dialog_compartment = new dialog_compartment.Dialog(this);

        this.reaction_page = new reactions.Page(document.getElementById("reaction-tab")!, this);
        this.metabolite_page = new metabolites.Page(document.getElementById("metabolite-tab")!, this);
        this.compartment_page = new compartments.Page(document.getElementById("compartment-tab")!, this);
        //this.stoichiometry_page = new stoichiometry.Page(document.getElementById("chemical-tab")!, this);
        this.settings_page = new settings.Page(document.getElementById("settings-tab")!, this);
        this.history_page = new history.Page(document.getElementById("history-tab")!, this);
        this.simulation_page = new simulation.Page(document.getElementById("simulation-tab")!, this);
    }
}

export function run(mm_cls: any, urls: any) {
    $.ajax({
        url: urls.get_reactions,
        context: document.body
    }).done(function(x: any) {
        app = new AppManager(mm_cls, x, urls);

        app.reaction_page.init();
        app.metabolite_page.init();
        app.compartment_page.init();
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

        $("#design-save").on("click", function() {
            app.dialog_save.show();
        });

        app.request_handler.endRequest($("#content"));
    });
}
