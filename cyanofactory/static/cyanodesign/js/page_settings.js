define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    let template = document.createElement('template');
    template.innerHTML = `
<div class="cyano-design-objective">
    <fieldset>
        <div class="form-group">
            <label class="control-label" for="cyano-design-objective-select">Main objective</label>
            <select class="cyano-design-objective-select" placeholder="Select a main objective"></select>
        </div>

       <div class="radio">
            <input type="radio" name="maxmingroup" class="radio-maximize" value="maximize-obj">
            <label for="radio-maximize">
                Maximize objective
            </label>
        </div>
        <div class="radio">
            <input type="radio" name="maxmingroup" class="radio-minimize" value="minimize-obj">
            <label for="radio-minimize">
                Minimize objective
            </label>
        </div>
    </fieldset>
</div>

<div class="panel panel-default">
    <div class="panel-heading">
        <h3 class="panel-title">Simulation Type</h3>
    </div>
    <div class="panel-body" class="simulation-type">
        <div class="radio">
            <input type="radio" name="simulationtype" class="radio-fba">
            <label for="radio-fba">
                Flux Balance Analysis
            </label>
        </div>
        <div class="radio">
            <input type="radio" name="simulationtype" class="radio-mba">
            <label for="radio-mba">
                Robustness Analysis
            </label>
        </div>
        <div class="radio">
            <input type="radio" name="simulationtype" class="radio-sa">
            <label for="radio-sa">
                Sensitivity Analysis
            </label>
        </div>
    </div>
</div>

<div class="panel panel-default" class="fba-settings-panel">
    <div class="panel-heading">
        <h3 class="panel-title">FBA settings</h3>
    </div>

    <div class="panel-body" class="settings-fba">
        <div>
            The main objective is always displayed in the simulation.
        </div>

        <div>
            <div class="radio">
                <input class="auto_flux" type="radio" name="flux" value="auto_flux">
                <label for="auto_flux">Display reactions with highest flux connected to the objective</label>
            </div>
            <div class="radio">
                <input class="manual_flux" type="radio" name="flux" value="manual_flux">
                <label for="manual_flux">Display selected reactions</label>
            </div>
        </div>

        <div class="cyano-design-reaction-filter">
            <fieldset>
                <div class="form-group">
                    <label class="control-label" for="design-objective-visible-combobox">Select visible reactions</label>
                    <input class="design-objective-visible-combobox">
                </div>
            </fieldset>
        </div>
    </div>
</div>

<div class="panel panel-default" class="sa-settings-panel">
    <div class="panel-heading">
        <h3 class="panel-title">Sensitivity/RA settings</h3>
    </div>

    <div class="panel-body" class="settings-mba">
        <fieldset>
            <div class="form-group">
                <label class="control-label" for="cyano-design-design-objective-select">Design objective</label>
                <select class="cyano-design-design-objective-select" placeholder="Select a design objective"></select>
            </div>
            <div class="form-group" class="sa-settings-target-function-panel">
                <label class="control-label" for="cyano-design-target-objective-select">Target reaction</label>
                <select class="cyano-design-target-objective-select" placeholder="Select a target reaction"></select>
            </div>
        </fieldset>
    </div>
</div>

<div class="panel panel-default" class="sa-settings-axis-panel">
    <div class="panel-heading">
        <h3 class="panel-title">Axis label</h3>
    </div>

    <div class="panel-body">
        <fieldset>
            <label class="control-label" for="mba-x-label">X-Axis</label>
            <input class="form-control" type="text" class="mba-x-label">
            <label class="control-label" for="mba-y-label">Y-Axis (Left)</label>
            <input class="form-control" type="text" class="mba-y-label">
            <label class="control-label" for="mba-y2-label">Y-Axis (Right)</label>
            <input class="form-control" type="text" class="mba-y2-label">
            Changes apply without running a new simulation
        </fieldset>
    </div>
</div>

`;
    class Page {
        constructor(where, app) {
            this.app = app;
            this.source_element = where;
            where.appendChild(template.content.cloneNode(true));
            this.main_obj_element = where.getElementsByClassName("cyano-design-objective-select")[0];
            this.design_obj_element = where.getElementsByClassName("cyano-design-design-objective-select")[0];
            this.target_obj_element = where.getElementsByClassName("cyano-design-target-objective-select")[0];
            this.maximize_element = where.getElementsByClassName("radio-maximize")[0];
            this.fba_sim_element = where.getElementsByClassName("radio-fba")[0];
            $(this.maximize_element).prop("checked", true);
            $(this.fba_sim_element).prop("checked", true);
            /*function updateSettingsVisibility() {
                $("#fba-settings-panel").hide();
                $("#sa-settings-panel").hide();
                $("#sa-settings-axis-panel").hide();
    
                if ($("#radio-fba").prop("checked")) {
                    $("#fba-settings-panel").show();
                } else if ($("#radio-mba").prop("checked")) {
                    $("#sa-settings-panel").show();
                    $("#sa-settings-target-function-panel").hide();
                    $("#sa-settings-axis-panel").show();
                } else if ($("#radio-sa").prop("checked")) {
                    $("#sa-settings-panel").show();
                    $("#sa-settings-target-function-panel").show();
                    $("#sa-settings-axis-panel").show();
                }
            }
    
            function setSimulationType(t) {
                if (t == "fba") {
                    $("#radio-fba").prop("checked", true);
                } else if (t == "mba") {
                    $("#radio-mba").prop("checked", true);
                } else if (t == "sa") {
                    $("#radio-sa").prop("checked", true);
                }
    
                updateSettingsVisibility();
            }
    
            $("#radio-fba").change(updateSettingsVisibility);
            $("#radio-mba").change(updateSettingsVisibility);
            $("#radio-sa").change(updateSettingsVisibility);*/
            /*For display list
                    design_objective_visible_combobox.selectize({
                        delimiter: '\x00',
                        labelField: "name",
                        valueField: "id",
                        searchField: ["name"],
                        plugins: ['drag_drop', 'remove_button', 'restore_on_backspace'],
                        maxItems: 20
                    });*/
            /*
            
                    $("#mba-x-label").change(function() {
                        updateLabels();
                    });
                     $("#mba-y-label").change(function() {
                        updateLabels();
                    });
                     $("#mba-y2-label").change(function() {
                        updateLabels();
                    });
            
             */
        }
        update() {
            const obj_options = {
                maxItems: 1,
                valueField: 'id',
                searchField: ['name', 'id'],
                options: this.app.model.reactions,
                render: {
                    item: function (item, escape) {
                        return "<div>" + escape(item.get_name_or_id()) + "</div>";
                    },
                    option: function (item, escape) {
                        return "<div>" + escape(item.get_name_or_id()) + "</div>";
                    }
                }
            };
            $(this.main_obj_element)["selectize"](obj_options);
            $(this.design_obj_element)["selectize"](obj_options);
            $(this.target_obj_element)["selectize"](obj_options);
        }
        getObjective() {
            return this.app.model.reaction.checked_get("id", this.main_obj_element["selectize"].getValue());
        }
        maximizeObjective() {
            return $(this.maximize_element).prop("checked");
        }
        getSimulationType() {
            if ($(this.fba_sim_element).prop("checked")) {
                return "fba";
            }
            else if ($(this.source_element.getElementsByClassName("radio-mba")[0]).prop("checked")) {
                return "mba";
            }
            else {
                return "ra";
            }
        }
        getFbaSettings() {
            return {
                highest_flux: true,
                selection: null
            };
        }
        getMbaSettings() {
            return {
                design_obj: null,
                target_reaction: null
            };
        }
        getAxisLabels() {
            return {
                x: "",
                y_left: "",
                y_right: ""
            };
        }
    }
    exports.Page = Page;
});
