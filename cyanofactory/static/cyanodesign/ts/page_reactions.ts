import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

let template = document.createElement('template');
template.innerHTML = `
<table class="cyano-reaction-list table table-striped table-hover">
    <thead>
        <tr>
            <th>Name</th>
            <th>Reaction</th>
            <th>Constraint</th>
            <th>Flux</th>
            <th>Active</th>
            <th></th>
        </tr>
    </thead>
</table>

<div class="cyano-reaction-buttons btn-group">
    <button type="button" class="btn btn-primary create-enzyme-button">
    Create new reaction
    </button>
    <button type="button" class="btn btn-primary dropdown-toggle" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
        <span class="caret"></span>
        <span class="sr-only">Toggle Dropdown</span>
    </button>
    <ul class="dropdown-menu">
        <li><a href="#" class="create-enzyme-bulk-button">Bulk add reactions</a></li>
    </ul>
</div>

<button type="button" class="btn btn-primary create-metabolite-button">
Create new metabolite
</button>
`;

let template_filter = document.createElement('template');
template_filter.innerHTML = `
<div class="col-sm-6">
<label for="cyano-list-filter">Filter reactions</label>
<select class="cyano-list-filter form-control combobox" multiple="multiple">
    <optgroup label="Enabled">
    <option selected="selected">Active</option>
    <option selected="selected">Inactive</option>
    </optgroup>
    <optgroup label="Constraints">
    <option selected="selected">Constraint</option>
    <option selected="selected">Unconstraint</option>
    </optgroup>
    <optgroup label="Reversible">
    <option selected="selected">Reversible</option>
    <option selected="selected">Irreversible</option>
    </optgroup>
</select>
</div>
<div class="col-sm-6">
    <div class="dataTables_filter">
    <div class="checkbox">
    <input class="cyano-regex" type="checkbox">
    <label for="cyano-regex">Search with RegExp</label>
    </div>
    </div>
</div>
`;

export class Page {
    readonly app: app.AppManager;
    readonly datatable: DataTables.Api;
    readonly source_element: HTMLElement;
    readonly table_element: HTMLElement;
    flux: any;

    constructor(where: HTMLElement, app: app.AppManager) {
        this.flux = {};

        this.source_element = where;
        where.appendChild(template.content.cloneNode(true));

        this.table_element = <HTMLElement>where.getElementsByClassName("cyano-reaction-list")[0]!;
        this.app = app;

        this.datatable = $(this.table_element).DataTable(<any>{
            "deferRender": true,
            columns: [
                    {},
                    {},
                    {},
                    {},
                    {},
                    {},
                    {}
                ],
            columnDefs: [
                {
                    "targets": 0,
                    "width": "25%",
                    "data": function(rowData: mm.Reaction) {
                        return rowData;
                    },
                    render: function(data: mm.Reaction) {
                        let e = $("<span>");

                        if (!data.enabled) {
                            e.addClass("cyano-enzyme-disabled");
                        }

                        return e.text(data.name).wrap("<p>").parent().html();
                    }
                },
                {
                    "targets": 1,
                    "width": "45%",
                    "orderable": false,
                    "data": function (rowData: mm.Reaction) {
                        return rowData;
                    },
                    render: function(data) {
                        return $(data.toHtml(app.model)).wrap("<p>").parent().html();
                    }
                },
                {
                    "targets": 2,
                    "width": "10%",
                    "orderable": false,
                    "searchable": false,
                    "data": function (rowData: mm.Reaction) {
                        return rowData.constraintsToString();
                    }
                },
                {
                    "targets": 3,
                    "width": "10%",
                    "orderable": false,
                    "searchable": false,
                    "data": function (rowData: mm.Reaction) {
                        if (rowData.id in self.flux) {
                            return self.flux[rowData.id];
                        }
                        return "";
                    }
                },
                {
                    "targets": 4,
                    "width": "12%",
                    "orderable": false,
                    "searchable": false,
                    "data": function (rowData: mm.Reaction) {
                        return rowData.enabled;
                    },
                    render: function(data: boolean, type, row, meta) {
                        if (type === 'copy') {
                            return data === true ? "Enabled" : "Disabled";
                        }

                        return "<div class='checkbox'> \
                        <input type='checkbox' id='enabled" + meta.row + "' " + (data ? "checked='checked'" : "") + "> \
                        <label for='enabled"  + meta.row + "'>Enabled</label> \
                        </div>";
                    }
                },
                {
                    "targets": 5,
                    "width": "8%",
                    "orderable": false,
                    "searchable": false,
                    "data": function (rowData: mm.Reaction, type, set, meta) {
                        return rowData;
                    },
                    "render": function(data, type, row, meta) {
                        return "<a class='btn btn-default btn-xs delete-button'>Delete</a>"
                    }
                },
                {
                    "targets": 6,
                    "visible": false,
                    "orderable": true,
                    "data": function (rowData, type, set, meta) {
                        return "No Pathway";
                        // FIXME return rowData.pathway.length == 0 ? "No Pathway" : rowData.pathway;
                    }
                }
            ],
            drawCallback: function ( settings ) {
                var api = this.api();
                var rows = api.rows( {page:'current'} ).nodes();
                var last=null;

                /*api.column(5, {page:'current'} ).data().each( function ( group, i ) {
                    if ( last !== group ) {
                        $(rows).eq( i ).before(
                            '<tr class="group"><td colspan="5">'+group+'</td></tr>'
                        );

                        last = group;
                    }
                } );*/
            },
            "displayLength": 25,
            "order": [[ 0, 'asc' ]],
            dom:
                "<'row'<'col-sm-6'l><'col-sm-6'f>>" +
                "<'row'>" +
                "<'row'<'col-sm-12'tr>>" +
                "<'row'<'col-sm-12'B>>" +
                "<'row'<'col-sm-6'i><'col-sm-6'p>>",
            buttons: [
                {
                    extend: 'copy',
                    exportOptions: {
                        columns: [0, 1, 2, 3],
                        modifier: {
                            selected: true
                        },
                        orthogonal: 'copy'
                    }
                },
                {
                    extend: 'csv',
                    exportOptions: {
                        columns: [0, 1, 2, 3],
                        modifier: {
                            selected: true
                        },
                        orthogonal: 'copy'
                    }
                }
            ]
        });

        /* Filter */
        where.children[0].children[1].appendChild(template_filter.content.cloneNode(true));

        const self: Page = this;
        (<any>$(where.getElementsByClassName("cyano-list-filter")[0])).multiselect({
            buttonClass: 'btn btn-default btn-xs',
            onChange: function(option, checked, select) {
                self.datatable.draw();
            },
            buttonText: function(options: HTMLOptionElement[], select) {
                if (options.length === 0) {
                    return 'No option selected';
                }
                else if (options.length === 6) {
                    return 'No filter applied';
                }
                else {
                    let labels = [];
                    let prev: HTMLOptionElement = null;
                    let omit_group: boolean = false;

                    let enabled: string[] = [];
                    for (let i = 0; i < 6; ++i) {
                        enabled.push("");
                    }

                    for (let i = 0; i < options.length; ++i) {
                        enabled[options[i].index] = options[i].text;
                    }

                    for (let i = 0; i < enabled.length; i+=2) {
                        let one_of: boolean = (enabled[i].length > 0 && enabled[i+1].length == 0) ||
                            (enabled[i].length == 0 && enabled[i+1].length > 0);
                        let none_of: boolean = (enabled[i].length == 0 && enabled[i+1].length == 0);

                        if (none_of) {
                            return "No option selected";
                        }

                        if (one_of) {
                            if (enabled[i].length > 0) {
                                labels.push(enabled[i]);
                            }
                            else if (enabled[i+1].length > 0) {
                                labels.push(enabled[i+1]);
                            }
                            labels.push("AND");
                        }
                    }

                    labels.pop();

                    return labels.join(' ') + '';
                }
            }
        });

        // Order by the grouping
        /*table_enzymes.delegate('tr.group', 'click', function() {
            var currentOrder = datatable_enzymes.order()[0];
            if ( currentOrder[0] === 2 && currentOrder[1] === 'asc' ) {
                datatable_enzymes.order( [ 5, 'desc' ] ).draw();
            }
            else {
                datatable_enzymes.order( [ 5, 'asc' ] ).draw();
            }
        } );*/

        $.fn.dataTable.ext.search.push(
            function( settings, data, dataIndex ) {
                if (settings.nTable == self.table_element) {
                    const arr = $(where).find(".cyano-list-filter").find("option").map(function () {
                        return (<any>this).selected;
                    }).get();

                    const d = self.datatable.data()[dataIndex];
                    const f = [
                        function(e) { return e.enabled },
                        function(e) { return !e.enabled },
                        function(e) { return e.isConstrained() },
                        function(e) { return !e.isConstrained() },
                        function(e) { return e.reversible },
                        function(e) { return !e.reversible }
                    ];

                    if (!(arr[0] || arr[1]) || !(arr[2] || arr[3]) || !(arr[4] || arr[5])) {
                        return false;
                    }

                    return ((arr[0] ? f[0](d) : !f[0](d)) || (arr[1] ? f[1](d) : !f[1](d))) &&
                            ((arr[2] ? f[2](d) : !f[2](d)) || (arr[3] ? f[3](d) : !f[3](d))) &&
                            ((arr[4] ? f[4](d) : !f[4](d)) || (arr[5] ? f[5](d) : !f[5](d)));
                }

                return true;
            }
        );

        where.getElementsByClassName("cyano-regex")[0].addEventListener("click", function() {
            self.datatable.search(self.datatable.search(), $(this).prop("checked"), true).draw();
        });

        /* Event handler */
        $(this.table_element).delegate('tr td:first-child', 'click', function() {
            // Reaction in 1st column was clicked
            let tr = $(this).closest("tr");
            if (tr.hasClass("group")) {
                return;
            }
            let row = self.datatable.row(tr);

            self.app.dialog_reaction.show(<mm.Reaction>row.data());
        });

        $(this.table_element).on("click", ".cyano-metabolite", function(event) {
            // Any Metabolite in 2nd column was clicked
            const met = self.app.model.metabolite.checked_get("id", this.dataset.id);
            self.app.dialog_metabolite.show(met);
        });

        // Enabled in 5th col
        $(this.table_element).delegate('tr td:nth-child(5) input', 'change', function() {
            let row = self.datatable.row($(this).closest("tr"));
            let reaction = (<mm.Reaction>row.data());

            reaction.enabled = ($(this).is(":checked"));

            app.history_manager.push({
                "type": "reaction",
                "op": "edit",
                "id": reaction.id,
                "object": {
                    "id": reaction.id,
                    "enabled": reaction.enabled
                }
            });

            self.invalidate(reaction);

            self.app.history_page.refresh();
        });

        // 5th col: delete button clicked
        $(this.table_element).on("click", ".delete-button", function() {
            let row = self.datatable.row($(this).closest("tr"));
            let reaction = (<mm.Reaction>row.data());

            self.app.dialog_reaction_delete.show(reaction);
        });
    }

    init() {
        this.datatable.clear();
        this.datatable.rows.add(this.app.model.reactions);

        for (let reaction of this.app.model.reactions) {
            reaction.updateMetaboliteReference(this.app.model);
        }

        this.refresh();
    }

    refresh() {
        this.datatable.sort();
        this.datatable.draw();
    }

    invalidate(reaction: mm.Reaction) {
        for (let met_id of reaction.getMetaboliteIds(this.app.model)) {
            let idx = this.app.model.metabolite.checked_index("id", met_id);
            this.app.metabolite_page.datatable.row(idx).invalidate("data");
        }

        this.app.metabolite_page.datatable.draw();

        (<any>this.datatable.row(this.app.model.reaction.checked_index("id", reaction.id))).invalidate();

        this.solve();
    }

    solve() {
        const solutions = [
            "Undefined",
            "Feasible",
            "Infeasible",
            "Not feasible",
            "Optimal",
            "Unbound"
        ];

        this.app.glpk_worker.onerror = (err) => {
            console.log(err);
        };

        this.app.glpk_worker.onmessage = (evt) => {
            //console.log(JSON.stringify(evt.data, null, 2));

            this.flux = {};
            const vars = evt.data.result.vars;
            for (const key in vars) {
                this.flux[key] = vars[key];
            }

            this.datatable.rows().invalidate();

            const fn = evt.data.result.status == 5 ? this.app.simulation_page.notifyInfo : this.app.simulation_page.notifyWarning;
            fn("The solution is " + solutions[evt.data.result.status - 1] + ". Flux of objective is " +
                evt.data.result.z.toFixed(4));
        };

        if (this.app.settings_page.getObjective() == "") {
            return;
        }

        let objective = {
            name: "z",
            direction: this.app.settings_page.maximizeObjective() ? 2 : 1,
            vars: [{
                name: this.app.settings_page.getObjective(),
                coef: 1.0
            }]
        };

        let subjectTo = [];
        let bounds = [];
        let transport_reacs: mm.Reaction[] = [];

        for (const met of this.app.model.metabolites) {
            subjectTo.push({
                name: met.id,
                vars: [
                ],
                bnds: { type: 5, ub: 0.0, lb: 0.0 }
            });

            if (met.isExternal(this.app.model)) {
                let r = new mm.Reaction();
                r.id = met.id + " <-> TRANSPORT";
                r.reversible = true;
                r.lower_bound = -10000.0;
                r.upper_bound = 10000.0;
                let mr = new mm.MetaboliteReference();
                mr.id = met.id;
                mr.stoichiometry = 1.0;
                r.products.push(mr);
                transport_reacs.push(r);
            }
        }

        for (const reac of this.app.model.reactions) {
            for (const s of reac.substrates) {
                subjectTo[this.app.model.metabolite.index("id", s.id)]["vars"].push({
                    name: reac.id,
                    coef: -s.stoichiometry
                });
            }

            for (const p of reac.products) {
                subjectTo[this.app.model.metabolite.index("id", p.id)]["vars"].push({
                    name: reac.id,
                    coef: p.stoichiometry
                });
            }

            if (reac.enabled) {
                bounds.push({
                    name: reac.id,
                    type: 4,
                    ub: reac.upper_bound,
                    lb: reac.lower_bound
                });
            } else {
                bounds.push({
                    name: reac.id,
                    type: 5,
                    ub: 0.0,
                    lb: 0.0
                });
            }
        }

        for (const reac of transport_reacs) {
            for (const p of reac.products) {
                subjectTo[this.app.model.metabolite.index("id", p.id)]["vars"].push({
                    name: reac.id,
                    coef: p.stoichiometry
                });
            }

            bounds.push({
                name: reac.id,
                type: 4,
                ub: reac.upper_bound,
                lb: reac.lower_bound
            });
        }

        this.app.glpk_worker.postMessage({
            name: 'LP',
            objective: objective,
            subjectTo: subjectTo,
            bounds: bounds
        });
    }
}
