const cyano_prefix = "../../cyano/node_modules/";

requirejs.config({
    paths: {
    }
});

// ignore datatables.net and jQuery, already included via script-tag
// causes havoc because jQuery is otherwise loaded twice
define('datatables.net', function() {
    return undefined;
});

define('jquery', function() {
    return $;
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
