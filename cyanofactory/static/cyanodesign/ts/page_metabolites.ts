import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

let template = document.createElement('template');
template.innerHTML = `
<table class="cyano-metabolite-list table table-striped table-hover">
    <thead>
        <tr>
            <th>Name</th>
            <th>Consumed by</th>
            <th>Produced by</th>
            <th>Is External</th>
        </tr>
    </thead>
</table>

<button type="button" class="btn btn-primary create-reaction-button">
Create new reaction
</button>
<button type="button" class="btn btn-primary create-metabolite-button">
Create new metabolite
</button>
<button type="button" class="btn btn-danger delete-metabolites-button">
Delete unused metabolites
</button>
`;

/*
<script id="filter_row_metabolites" type="text/plain">
<div class="col-sm-6">
<label for="cyano-metabolite-list-filter-reactions">Filter metabolites</label>
<select id="cyano-metabolite-list-filter-reactions" class="form-control combobox" multiple="multiple">
    <option selected="selected">Is internal</option>
    <option selected="selected">Is external</option>
</select>
</div>
<div class="col-sm-6">
    <div class="dataTables_filter">
    <div class="checkbox">
    <input id="metabolite_regex" type="checkbox">
    <label for="metabolite_regex">Search with RegExp</label>
    </div>
    </div>
</div>
</script>
 */

/*

    var metabolite_filter = [
        function(e) { return !e.isExternal() },
        function(e) { return e.isExternal() }
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

        this.table_element = <HTMLElement>where.getElementsByClassName("cyano-metabolite-list")[0]!;
        this.app = app;

        this.datatable = $(this.table_element).DataTable(<any>{
            "deferRender": true,
            columns: [
                    {},
                    {},
                    {},
                    {}
                ],
            columnDefs: [
                {
                    "targets": 0,
                    "data": function (rowData, type, set, meta) {
                        return rowData.name;
                    }
                },
                {
                    "targets": 1,
                    "orderable": false,
                    "data": function (rowData, type, set, meta) {
                        return rowData.consumed;
                    },
                    "render": function(data, type, row, meta ) {
                        return data.map(function(arg) {
                            var e = $("<span></span>").addClass("cyano-enzyme");
                            if (!arg.enabled) {
                                e.addClass("cyano-enzyme-disabled");
                            }
                            return e.append(arg.name).wrap("<p>").parent().html();
                        }).join(", ");
                    }
                },
                {
                    "targets": 2,
                    "orderable": false,
                    "data": function (rowData, type, set, meta) {
                        return rowData.produced;
                    },
                    "render": function(data, type, row, meta ) {
                        return data.map(function(arg) {
                            var e = $("<span></span>").addClass("cyano-enzyme");
                            if (!arg.enabled) {
                                e.addClass("cyano-enzyme-disabled");
                            }
                            return e.append(arg.name).wrap("<p>").parent().html();
                        }).join(", ");
                    }
                },
                {
                    "targets": 3,
                    "orderable": false,
                    "searchable": false,
                    "data": function (rowData, type, set, meta) {
                        return rowData.isExternal();
                    },
                    render: function(data, type, row, meta) {
                        if (type === 'copy') {
                            return data === true ? "External" : "Internal";
                        }

                        return "<div class='checkbox'> \
                        <input type='checkbox' id='external" + meta.row + "' " + (data ? "checked='checked'" : "") + "> \
                        <label for='external"  + meta.row + "'>External</label> \
                        </div>";
                    }
                }
            ],
            "displayLength": 25,
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

        /* Event handler */
        const self: Page = this;

        // Tooltip on reaction hover
        $(this.table_element).on({
            mouseenter: function () {
                $(this).data("toggle", "tooltip");
                $(this).data("placement", "top");

                let reaction: mm.Reaction = self.app.model.reaction.checked_get("name", $(this).text());
                let text: string = reaction.toString();

                $(this).prop("title", text);
                $(this)["tooltip"]("show");
            },
            mouseleave: function () {
            }
        }, ".cyano-enzyme");

        // 1st column clicked
        $(this.table_element).delegate('tr td:first-child', 'click', function() {
            let row = self.datatable.row($(this).closest("tr"));
            app.dialog_metabolite.show(<mm.Metabolite>row.data());
        });

        /*
// 4th col
table_metabolites.delegate('tr td:nth-child(4) input', 'change', function() {
    var row = datatable_metabolites.row($(this).closest("tr"));

    row.data().external = $(this).is(":checked");

    command_list.push({
        "type": "metabolite",
        "op": "edit",
        "id": row.data().id,
        "object": {
            "id": row.data().id,
            "name": row.data().name,
            "external": row.data().external
        }
    });

    row.data().invalidate();
});

// Enzyme in 2nd or 3rd col
table_metabolites.on("click", ".cyano-enzyme", function(event) {
    var enzyme = Enzyme.indexByName($(this).text());
    if (enzyme >= 0) {
        showEditEnzymeDialog(model.reactions[enzyme], false);
    }
});


        $(".create-enzyme-button").click(function (event) {
            showAddEnzymeDialog();
        });
        $(".create-metabolite-button").click(function (event) {
            showAddMetaboliteDialog();
        });

        $(".delete-metabolites-button").click(function (event) {
            $("#dialog-delete-metabolites").modal('show');
        });

        $("#dialog-delete-metabolites").find(".btn-primary").click(function() {
            model.metabolites.filter(function(m) {
                return m.isUnused();
            }).forEach(function(m) {
                m.remove();
                m.removeFromList();

                command_list.push({
                    "type": "metabolite",
                    "op": "delete",
                    "id": m.id,
                    "object": {}
                });
            });
            datatable_metabolites.draw();

            $("#dialog-delete-metabolites").modal("hide");
        });

*/
    }

    update() {
        this.datatable.clear();
        this.datatable.rows.add(this.app.model.metabolites);
        this.datatable.draw();
    }
}
