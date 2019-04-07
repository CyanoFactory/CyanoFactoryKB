const cyano_prefix = "../../cyano/node_modules/";

requirejs.config({
    paths: {
        jquery: cyano_prefix + 'jquery/dist/jquery',
        "datatables.net": cyano_prefix + 'datatables.net/js/jquery.dataTables'
    }
});

define([
    "./app",
    "./metabolic_model",
    "./page_reactions",
    "./page_metabolites",
    "./page_stoichiometry",
    "./page_settings",
    "./page_simulation",
    "./dialog_reaction",
    "./dialog_reaction_bulkadd",
    "./dialog_metabolite"
], function(app, mm) {
    MetabolicModel = mm;
    startup(app.run);
});
