define(["require", "exports", "./dialog_helper", "jquery", "datatables.net", "selectize"], function (require, exports, dialog_helper_1, $) {
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

<div class="panel panel-default fba-settings-panel">
    <div class="panel-heading">
        <h3 class="panel-title">FBA settings</h3>
    </div>

    <div class="panel-body settings-fba">
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

<div class="panel panel-default sa-settings-panel">
    <div class="panel-heading">
        <h3 class="panel-title">Sensitivity/RA settings</h3>
    </div>

    <div class="panel-body settings-mba">
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

<div class="panel panel-default sa-settings-axis-panel">
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

<div class="panel panel-default other-settings-axis-panel">
    <div class="panel-heading">
        <h3 class="panel-title">Other</h3>
    </div>

    <div class="panel-body">
        <fieldset>
            <label class="control-label" for="exchange-reaction-label">Add exchange reactions:</label>
            <select class="exchange-reaction-setting">
                <option value="auto" selected>Automatic (Recommended)</option>
                <option value="yes">Yes</option>
                <option value="no">No</option>
            </select>
            When the model does not provide exchange reactions for external metabolites, select "Yes", otherwise "No". Keep it on "Auto" if you are not sure.
        </fieldset>
    </div>
</div>

`;
    class Page {
        constructor(where, app) {
            this.app = app;
            this.source_element = where;
            where.appendChild(template.content.cloneNode(true));
            this.main_obj = dialog_helper_1.ElementWrapper.byClass("cyano-design-objective-select", this.source_element);
            this.design_obj = dialog_helper_1.ElementWrapper.byClass("cyano-design-design-objective-select", this.source_element);
            this.target_obj = dialog_helper_1.ElementWrapper.byClass("cyano-design-target-objective-select", this.source_element);
            this.exchange_reaction = dialog_helper_1.ElementWrapper.byClass("exchange-reaction-setting", this.source_element);
            this.maximize = dialog_helper_1.ElementWrapper.byClass("radio-maximize", this.source_element);
            this.fba_sim = dialog_helper_1.ElementWrapper.byClass("radio-fba", this.source_element);
            this.mba_sim = dialog_helper_1.ElementWrapper.byClass("radio-mba", this.source_element);
            this.sa_sim = dialog_helper_1.ElementWrapper.byClass("radio-sa", this.source_element);
            this.maximize.value = true;
            this.fba_sim.value = true;
            $(this.fba_sim.element).change(() => Page.updateSettingsVisibility(this));
            $(this.mba_sim.element).change(() => Page.updateSettingsVisibility(this));
            $(this.sa_sim.element).change(() => Page.updateSettingsVisibility(this));
            Page.updateSettingsVisibility(this);
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
        init() {
            const obj_options = {
                maxItems: 1,
                valueField: 'id',
                searchField: ['name', 'id'],
                render: {
                    item: function (item, escape) {
                        return "<div>" + escape(item.get_name_or_id()) + "</div>";
                    },
                    option: function (item, escape) {
                        return "<div>" + escape(item.get_name_or_id()) + "</div>";
                    }
                }
            };
            $(this.main_obj.element).selectize(obj_options);
            $(this.design_obj.element).selectize(obj_options);
            $(this.target_obj.element).selectize(obj_options);
            $(this.exchange_reaction.element).selectize();
            const self = this;
            let main_obj_selectize = this.main_obj.element.selectize;
            main_obj_selectize.on('change', function () {
                self.app.simulation_page.solve();
            });
            let exchange_reaction = this.exchange_reaction.element.selectize;
            exchange_reaction.on('change', function () {
                self.app.simulation_page.solve();
            });
            this.refresh();
        }
        refresh(reactions = null) {
            let main_obj_selectize = this.main_obj.element.selectize;
            let design_obj_selectize = this.design_obj.element.selectize;
            let target_obj_selectize = this.target_obj.element.selectize;
            if (reactions == null) {
                const fn_c = (s) => {
                    s.clearOptions();
                };
                const fn = (r, s) => {
                    s.addOption(r);
                };
                fn_c(main_obj_selectize);
                fn_c(design_obj_selectize);
                fn_c(target_obj_selectize);
                for (const reac of this.app.model.reactions) {
                    fn(reac, main_obj_selectize);
                    fn(reac, design_obj_selectize);
                    fn(reac, target_obj_selectize);
                }
            }
            else {
                const fn = (r, s) => {
                    if (r.id == s.getValue()) {
                        s.clear();
                    }
                    s.removeOption(r.id);
                };
                for (const reac of reactions) {
                    fn(reac, main_obj_selectize);
                    fn(reac, design_obj_selectize);
                    fn(reac, target_obj_selectize);
                }
            }
        }
        getObjective() {
            return this.main_obj.element.selectize.getValue();
        }
        getDesignObjective() {
            return this.design_obj.element.selectize.getValue();
        }
        getTargetObjective() {
            return this.target_obj.element.selectize.getValue();
        }
        getCreateExchangeReactions() {
            const val = this.exchange_reaction.element.selectize.getValue();
            if (val == "yes") {
                return true;
            }
            else if (val == "no") {
                return false;
            }
            // Autodetect
            // When the model follows the cobra convention the exchange reacs are prefixed with "EX_"
            for (const met of this.app.model.getExternalMetabolites()) {
                if (met.id.startsWith("M_")) {
                    if (this.app.model.reaction.has("id", "EX_" + met.id.substr(2))) {
                        return false;
                    }
                    if (this.app.model.reaction.has("id", "R_EX_" + met.id.substr(2))) {
                        return false;
                    }
                }
                if (this.app.model.reaction.has("id", "EX_" + met.id)) {
                    return false;
                }
                if (this.app.model.reaction.has("id", "R_EX_" + met.id)) {
                    return false;
                }
            }
            return true;
        }
        maximizeObjective() {
            return this.maximize.value;
        }
        getSimulationType() {
            if (this.fba_sim.value) {
                return "fba";
            }
            else if (this.mba_sim.value) {
                return "mba";
            }
            else if (this.sa_sim.value) {
                return "sa";
            }
            // unreachable path
            throw new Error("getSimulationType: BUG! No type selected");
        }
        setSimulationType(t) {
            if (t == "fba") {
                this.fba_sim.value = true;
            }
            else if (t == "mba") {
                this.mba_sim.value = true;
            }
            else if (t == "sa") {
                this.sa_sim.value = true;
            }
            Page.updateSettingsVisibility(this);
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
        static updateSettingsVisibility(self) {
            $(self.source_element).find(".fba-settings-panel").hide();
            $(self.source_element).find(".sa-settings-panel").hide();
            $(self.source_element).find(".sa-settings-axis-panel").hide();
            if (self.fba_sim.value) {
                $(self.source_element).find(".fba-settings-panel").show();
            }
            else if (self.mba_sim.value) {
                $(self.source_element).find(".sa-settings-panel").show();
                $(self.source_element).find(".sa-settings-target-function-panel").hide();
                $(self.source_element).find(".sa-settings-axis-panel").show();
            }
            else if (self.sa_sim.value) {
                $(self.source_element).find(".sa-settings-panel").show();
                $(self.source_element).find(".sa-settings-target-function-panel").show();
                $(self.source_element).find(".sa-settings-axis-panel").show();
            }
        }
        ;
    }
    exports.Page = Page;
});
