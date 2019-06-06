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

// already loaded via scripts tag
define('selectize', function() {
    return undefined;
});

define([
    "./app",
    "./metabolic_model",
    "./urls",
    "./page_reactions",
    "./page_metabolites",
    "./page_stoichiometry",
    "./page_settings",
    "./page_simulation",
    "./page_history",
    "./dialog_reaction",
    "./dialog_reaction_bulkadd",
    "./dialog_save",
    "./dialog_metabolite",
    "./request_handler",
], function(app, mm, url) {
    startup(app.run, mm, url);
});
