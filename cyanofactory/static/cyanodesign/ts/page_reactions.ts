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
        <li><a href="#" id="create-enzyme-bulk-button">Bulk add reactions</a></li>
    </ul>
</div>

<button type="button" class="btn btn-primary create-metabolite-button">
Create new metabolite
</button>
`;

/*

<script id="filter_row" type="text/plain">
<div class="col-sm-6">
<label for="cyano-enzyme-list-filter-reactions">Filter reactions</label>
<select id="cyano-enzyme-list-filter-reactions" class="form-control combobox" multiple="multiple">
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
    <input id="enzyme_regex" type="checkbox">
    <label for="enzyme_regex">Search with RegExp</label>
    </div>
    </div>
</div>
</script>


 */

/*

    var reaction_filter = [
        function(e) { return e.enabled },
        function(e) { return !e.enabled },
        function(e) { return e.isConstrained() },
        function(e) { return !e.isConstrained() },
        function(e) { return e.reversible },
        function(e) { return !e.reversible }
    ];
 */

export class Page {
    readonly app: app.AppManager;
    readonly datatable: DataTables.Api;
    readonly source_element: HTMLElement;
    readonly table_element: HTMLElement;

    constructor(where: HTMLElement, app: app.AppManager) {
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
                        return $(data.toHTML(app.model)).wrap("<p>").parent().html();
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
                    "targets": 4,
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
                    "targets": 5,
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

        /* Register event handler */
        const self: Page = this;

        $(this.table_element).delegate('tr td:first-child', 'click', function() {
            // Reaction in 1st column was clicked
            let tr = $(this).closest("tr");
            if (tr.hasClass("group")) {
                return;
            }
            let row = self.datatable.row(tr);

            self.app.dialog_reaction.show();//(row.data(), false);
        });


        /*table_enzymes.on("click", ".cyano-metabolite", function(event) {
            // Any Metabolite in 2nd column was clicked
            var met = Metabolite.indexByName($(this).text());
            if (met >= 0) {
                showEditMetaboliteDialog(model.metabolites[met]);
            }
        });*/


        // First col
        /*table_enzymes.delegate('tr td:first-child', 'click', function() {
            var tr = $(this).closest("tr");
            if (tr.hasClass("group")) {
                return;
            }
            var row = datatable_enzymes.row(tr);

            showEditEnzymeDialog(row.data(), false);
        });



        // Enabled in 4th col
        table_enzymes.delegate('tr td:nth-child(4) input', 'change', function() {
            var row = datatable_enzymes.row($(this).closest("tr"));

            row.data().enabled = ($(this).is(":checked"));

            command_list.push({
                "type": "reaction",
                "op": "edit",
                "id": row.data().id,
                "object": {
                    "id": row.data().id,
                    "enabled": row.data().enabled
                }
            });

            row.data().invalidate();
        });*/

    }

    update() {
        this.datatable.clear();
        this.datatable.rows.add(this.app.model.reactions);
        this.datatable.draw();

        for (let reaction of this.app.model.reactions) {
            reaction.updateMetaboliteReference(this.app.model);
        }

        /*var old_design_objective = design_objective;
        model.reactions.forEach(function(item) {
            cyano_design_objective_select[0].selectize.addOption(item);
            cyano_design_design_objective_select[0].selectize.addOption(item);
            design_objective_visible_combobox[0].selectize.addOption(item);
            cyano_design_target_objective_select[0].selectize.addOption(item);
        });
        cyano_design_objective_select[0].selectize.refreshOptions();
        cyano_design_design_objective_select[0].selectize.refreshOptions();
        design_objective_visible_combobox[0].selectize.refreshOptions();
        cyano_design_target_objective_select[0].selectize.refreshOptions();

        if (old_design_objective !== undefined) {
            cyano_design_objective_select[0].selectize.setValue(old_design_objective.name, false);
        }*/
    }
}
