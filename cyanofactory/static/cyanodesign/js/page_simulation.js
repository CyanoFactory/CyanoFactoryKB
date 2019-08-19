define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    let template = document.createElement('template');
    template.innerHTML = `
<div class="checkbox">
    <input type="checkbox" name="remember-simulation" id="remember-simulation">
    <label for="remember-simulation">Combine results with previous simulation</label>
</div>
<button type="button" class="design-submit btn btn-primary">Run simulation</button>
<div class="export-button btn-group">
    <button type="button" class="btn btn-primary dropdown-toggle" data-toggle="dropdown" aria-expanded="false">
    Export <span class="caret"></span>
    </button>
    <ul class="dropdown-menu" role="menu">
    <!--<li><a id="export-csv" href="#">As flux list</a></li>-->
    <li><a id="export-png" href="#">As image</a></li>
    <li><a id="export-svg" href="#">As vector graphic</a></li>
    </ul>
</div>

<div class="simulation-result">
</div>

<div class="visual-graph">
</div>

<div class="visual-fba">
</div>

<table class="cyano-flux-list table table-striped table-hover">
    <thead>
        <tr>
            <th>Name</th>
            <th>Flux</th>
        </tr>
    </thead>
</table>
`;
    /*
    <script id="filter_row_flux" type="text/plain">
    <div class="col-sm-3">
        <div class="form-group">
            <label class="control-label" for="filter-flux-min">Min Flux:</label>
            <input id="filter-flux-min" class="form-control" type="text">
        </div>
    </div>
    <div class="col-sm-3">
            <label class="control-label" for="filter-flux-max">Max Flux:</label>
            <input id="filter-flux-max" class="form-control" type="text">
    </div>
    <!--
    <label for="cyano-metabolite-list-filter-reactions">Filter metabolites</label>
    <select id="cyano-metabolite-list-filter-reactions" class="form-control combobox" multiple="multiple">
        <option selected="selected">Is internal</option>
        <option selected="selected">Is external</option>
    </select>-->
    </div>
    <div class="col-sm-6">
        <div class="dataTables_filter">
        <div class="checkbox">
        <input id="flux_regex" type="checkbox">
        <label for="flux_regex">Search with RegExp</label>
        </div>
        </div>
    </div>
    </script>
     */
    class Page {
        constructor(where, app) {
            this.last_sim_type = "";
            this.last_sim_flux = 0;
            this.simulation_chart = null;
            this.is_dragging = false;
            this.source_element = where;
            where.appendChild(template.content.cloneNode(true));
            this.table_element_flux = where.getElementsByClassName("cyano-flux-list")[0];
            this.simulation_result_element = where.getElementsByClassName("simulation-result")[0];
            this.visual_graph_element = where.getElementsByClassName("visual-graph")[0];
            this.visual_fba_element = where.getElementsByClassName("visual-fba")[0];
            this.app = app;
            this.datatable_flux = $(this.table_element_flux).DataTable({
                "deferRender": true,
                "displayLength": 25,
                "order": [[1, 'desc']],
                "language": {
                    "emptyTable": "The fluxes are displayed here after running a simulation"
                },
                columnDefs: [
                    {
                        "targets": 0
                    },
                    {
                        "targets": 1
                    }
                ],
                dom: "<'row'<'col-sm-6'l><'col-sm-6'f>>" +
                    "<'row'>" +
                    "<'row'<'col-sm-12'tr>>" +
                    "<'row'<'col-sm-12'B>>" +
                    "<'row'<'col-sm-5'i><'col-sm-7'p>>",
                buttons: [
                    {
                        extend: 'copy',
                        exportOptions: {
                            modifier: {
                                selected: true
                            },
                            orthogonal: 'copy'
                        }
                    },
                    {
                        extend: 'csv',
                        exportOptions: {
                            modifier: {
                                selected: true
                            },
                            orthogonal: 'copy'
                        }
                    }
                ]
            });
            let self = this;
            $(this.source_element).find(".design-submit").click(function (event) {
                self.simulate();
            });
            /*FIXME$("#visual_fba").on("click", "g.node", function() {
                var idx = Metabolite.indexByName($(this).children("text").text());
                if (!isDragging & idx > -1) {
                    showEditMetaboliteDialog(model.metabolites[idx]);
                }
            });
    
            $("#visual_fba").on("click", "g.edge text", function() {
                // remove flux value
                var idx = Enzyme.indexByName($(this).text().replace(/ \(-?[0-9]+(\.[0-9]+)?\)$/, ""));
                if (!isDragging & idx > -1) {
                    showEditEnzymeDialog(model.reactions[idx], false);
                }
            });
    
            $("#visual_fba").on("mousedown", "g", function(event) {
                isDragging = false;
                $(this).data('page', {x: event.pageX, y: event.pageY})
            });
            $("#visual_fba").on("mousemove", "g", function(event) {
                var p = $(this).data('page');
                if (p !== undefined) {
                    if (Math.abs(p.x - event.pageX) > 4 ||
                        Math.abs(p.y - event.pageY) > 4) {
                        isDragging = true;
                    }
                }
            });
    
            $("#filter-flux-min").change(function() { datatable_flux.draw(); });
            $("#filter-flux-max").change(function() { datatable_flux.draw(); });
    */
        }
        init() {
        }
        updateLabels() {
            if (this.simulation_chart == null) {
                return;
            }
            this.simulation_chart.axis.labels({
                x: $("#mba-x-label").val(),
                y: $("#mba-y-label").val(),
                y2: $("#mba-y2-label").val()
            });
        }
        ;
        notifyInfo(text) {
            document.getElementById("wedesign-notify-box").innerHTML = '<div class="alert alert-info" role="alert">\
            <span class="glyphicon glyphicon-info-sign" aria-hidden="true"></span>\
            <span class="sr-only">Info:</span>' + text + '</div>';
        }
        notifyWarning(text) {
            document.getElementById("wedesign-notify-box").innerHTML = '<div class="alert alert-warning" role="alert">\
            <span class="glyphicon glyphicon-warning-sign" aria-hidden="true"></span>\
            <span class="sr-only">Info:</span>' + text + '</div>';
        }
        notifyError(text) {
            document.getElementById("wedesign-notify-box").innerHTML = '<div class="alert alert-danger" role="alert">\
            <span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>\
            <span class="sr-only">Error:</span>' + text + '</div>';
        }
        simulate() {
            let symtype = this.app.settings_page.getSimulationType();
            this.last_sim_type = symtype;
            if (symtype == "fba") {
                let graph = Viz(this.createGraph(this.app.reaction_page.flux), "svg", "dot");
                /*if (simulation_result["solution"] == "Optimal") {
                    var obj = this.app.settings_page.getObjective();
                    this.notifyInfo("The solution is " + simulation_result["solution"] + ". Flux of objective is " + simulation_result["flux"][obj].toFixed(4));
                } else {
                    this.notifyWarning("The solution is " + simulation_result["solution"] + ". Check if your constraints are too strict.");
                }
    */
                $(this.visual_graph_element).hide();
                $(this.visual_fba_element).show();
                this.visual_fba_element.innerHTML = graph;
                $(this.visual_fba_element).attr("width", "100%").attr("height", "400px");
                let svgPan = svgPanZoom('.visual-fba > svg', { minZoom: 0.1, fit: false });
                svgPan.zoom(1);
                this.datatable_flux.clear();
                for (const reac of this.app.model.reactions) {
                    this.datatable_flux.row.add([reac.get_name_or_id(), this.app.reaction_page.flux[reac.id]]);
                }
                this.datatable_flux.draw();
                /*$("svg").find(".edge text").each(function() {
                        if ($(this).text().indexOf(cyano_design_objective_select.val() + " ") == 0) {
                            var new_x = ($("svg").width() / 2) - $(this).attr("x");
                            var new_y = ($("svg").height() / 2) - $(this).attr("y");
                            svgPan.pan({x: new_x, y: new_y});
                        }
                    });*/
            }
            else if (symtype == "mba") {
                $(this.visual_graph_element).show();
                $(this.visual_fba_element).hide();
                if (!$("#remember-simulation").prop("checked") || this.simulation_chart === undefined) {
                    var chart = {
                        bindto: ".visual-graph",
                        data: {
                            x: 'x',
                            columns: simulation_result["graph"],
                            type: 'bar',
                            axes: {} /*,
                            onclick: chart_clicked*/
                        },
                        axis: {
                            x: {
                                label: {
                                    position: "outer-center"
                                },
                                type: 'category'
                            },
                            y: {
                                label: {
                                    position: "outer-middle"
                                }
                            }
                        }
                    };
                    chart["data"]["axes"][this.app.settings_page.getObjective()] = "y";
                    this.simulation_chart = c3.generate(chart);
                }
                else {
                    this.simulation_chart.load({
                        columns: [
                        //FIXME[cyano_design_design_objective_select.val()].concat(simulation_result["graph"][1])
                        ],
                        type: 'bar'
                    });
                }
                this.last_sim_flux = simulation_result["flux"];
                this.updateLabels();
            }
            else if (symtype == "sa") {
                $(this.visual_graph_element).show();
                $(this.visual_fba_element).hide();
                let chart = {
                    bindto: '.visual-graph',
                    data: {
                        x: 'x',
                        columns: simulation_result["graph"],
                        type: 'bar',
                        axes: {} /*,
                        onclick: chart_clicked*/
                    },
                    axis: {
                        x: {
                            label: {
                                position: "outer-center"
                            },
                            type: 'category'
                        },
                        y: {
                            label: {
                                position: "outer-middle"
                            }
                        },
                        y2: {
                            label: {
                                position: "outer-middle"
                            },
                            show: true
                        }
                    }
                };
                chart["data"]["axes"][this.app.settings_page.getObjective()] = "y";
                //FIXMEchart["data"]["axes"][cyano_design_design_objective_select.val()] = "y2";
                this.simulation_chart = c3.generate(chart);
                this.last_sim_flux = simulation_result["flux"];
                this.updateLabels();
            }
        }
        createGraph(flux) {
            let graph = 'strict digraph "" {\n' +
                'graph [overlap=False,\n' +
                'rankdir=LR,\n' +
                'splines=True\n' +
                '];\n' +
                'node [colorscheme=pastel19,\n' +
                'label="\\N",' +
                'style=filled\n' +
                '];\n';
            const reac_offset = this.app.model.metabolites.length;
            // Create sparse array
            let edges = [];
            // A_ext -> A
            let i = 0;
            for (const met of this.app.model.metabolites) {
                graph += "" + i + "[color=" + ((i % 10) + 1) + ",\n" +
                    'label=' + met.name + ',\n' +
                    'shape=oval];\n';
                ++i;
            }
            // Pen scaling factor calculation
            let max_flux = Number.MIN_SAFE_INTEGER;
            for (const f in flux) {
                max_flux = Math.max(flux[f] < 0.0 ? -flux[f] : flux[f], max_flux);
            }
            max_flux = 10 / max_flux;
            const obj = this.app.settings_page.getObjective();
            for (const reac of this.app.model.reactions) {
                for (const p of reac.products) {
                    for (const s of reac.substrates) {
                        let f = flux[reac.id];
                        edges.push({
                            label: reac.name,
                            flux: f,
                            color: f < 0.0 ? "red" : f > 0.0 ? "green" : "black",
                            penwidth: f == 0.0 ? 1 : f * max_flux,
                            left: this.app.model.metabolite.index("id", s.id),
                            right: this.app.model.metabolite.index("id", p.id),
                            reverse: reac.reversible,
                            obj: reac.id == obj
                        });
                    }
                }
            }
            i = reac_offset;
            for (const edge of edges) {
                graph += edge.left + " -> " + edge.right + " [color=" + edge.color + ",\n" +
                    'label="' + edge.label + ' (' + edge.flux + ')",\n' +
                    'penwidth=' + edge.penwidth + (edge.obj ? ',style=dashed' : '') + '];\n';
                if (edge.reverse) {
                    graph += edge.right + " -> " + edge.left + " [\n" +
                        'label="' + edge.label + '",\n' +
                        '];\n';
                }
                ++i;
            }
            graph += '}\n';
            return graph;
        }
    }
    exports.Page = Page;
});
