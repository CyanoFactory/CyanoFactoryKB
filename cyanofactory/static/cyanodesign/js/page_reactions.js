define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
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
    class Page {
        constructor(where, app) {
            this.source_element = where;
            where.appendChild(template.content.cloneNode(true));
            this.table_element = where.getElementsByClassName("cyano-reaction-list")[0];
            this.app = app;
            this.datatable = $(this.table_element).DataTable({
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
                        "data": function (rowData) {
                            return rowData;
                        },
                        render: function (data) {
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
                        "data": function (rowData) {
                            return rowData;
                        },
                        render: function (data) {
                            return $(data.toHtml(app.model)).wrap("<p>").parent().html();
                        }
                    },
                    {
                        "targets": 2,
                        "width": "10%",
                        "orderable": false,
                        "searchable": false,
                        "data": function (rowData) {
                            return rowData.constraintsToString();
                        }
                    },
                    {
                        "targets": 3,
                        "width": "12%",
                        "orderable": false,
                        "searchable": false,
                        "data": function (rowData) {
                            return rowData.enabled;
                        },
                        render: function (data, type, row, meta) {
                            if (type === 'copy') {
                                return data === true ? "Enabled" : "Disabled";
                            }
                            return "<div class='checkbox'> \
                        <input type='checkbox' id='enabled" + meta.row + "' " + (data ? "checked='checked'" : "") + "> \
                        <label for='enabled" + meta.row + "'>Enabled</label> \
                        </div>";
                        }
                    },
                    {
                        "targets": 4,
                        "width": "8%",
                        "orderable": false,
                        "searchable": false,
                        "data": function (rowData, type, set, meta) {
                            return rowData;
                        },
                        "render": function (data, type, row, meta) {
                            return "<a class='btn btn-default btn-xs delete-button'>Delete</a>";
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
                drawCallback: function (settings) {
                    var api = this.api();
                    var rows = api.rows({ page: 'current' }).nodes();
                    var last = null;
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
                "order": [[0, 'asc']],
                dom: "<'row'<'col-sm-6'l><'col-sm-6'f>>" +
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
            const self = this;
            $(where.getElementsByClassName("cyano-list-filter")[0]).multiselect({
                buttonClass: 'btn btn-default btn-xs',
                onChange: function (option, checked, select) {
                    self.datatable.draw();
                },
                buttonText: function (options, select) {
                    if (options.length === 0) {
                        return 'No option selected';
                    }
                    else if (options.length === 6) {
                        return 'No filter applied';
                    }
                    else {
                        let labels = [];
                        let prev = null;
                        let omit_group = false;
                        let enabled = [];
                        for (let i = 0; i < 6; ++i) {
                            enabled.push("");
                        }
                        for (let i = 0; i < options.length; ++i) {
                            enabled[options[i].index] = options[i].text;
                        }
                        for (let i = 0; i < enabled.length; i += 2) {
                            let one_of = (enabled[i].length > 0 && enabled[i + 1].length == 0) ||
                                (enabled[i].length == 0 && enabled[i + 1].length > 0);
                            let none_of = (enabled[i].length == 0 && enabled[i + 1].length == 0);
                            if (none_of) {
                                return "No option selected";
                            }
                            if (one_of) {
                                if (enabled[i].length > 0) {
                                    labels.push(enabled[i]);
                                }
                                else if (enabled[i + 1].length > 0) {
                                    labels.push(enabled[i + 1]);
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
            $.fn.dataTable.ext.search.push(function (settings, data, dataIndex) {
                if (settings.nTable == self.table_element) {
                    const arr = $(where).find(".cyano-list-filter").find("option").map(function () {
                        return this.selected;
                    }).get();
                    const d = self.datatable.data()[dataIndex];
                    const f = [
                        function (e) { return e.enabled; },
                        function (e) { return !e.enabled; },
                        function (e) { return e.isConstrained(); },
                        function (e) { return !e.isConstrained(); },
                        function (e) { return e.reversible; },
                        function (e) { return !e.reversible; }
                    ];
                    if (!(arr[0] || arr[1]) || !(arr[2] || arr[3]) || !(arr[4] || arr[5])) {
                        return false;
                    }
                    return ((arr[0] ? f[0](d) : !f[0](d)) || (arr[1] ? f[1](d) : !f[1](d))) &&
                        ((arr[2] ? f[2](d) : !f[2](d)) || (arr[3] ? f[3](d) : !f[3](d))) &&
                        ((arr[4] ? f[4](d) : !f[4](d)) || (arr[5] ? f[5](d) : !f[5](d)));
                }
                return true;
            });
            where.getElementsByClassName("cyano-regex")[0].addEventListener("click", function () {
                self.datatable.search(self.datatable.search(), $(this).prop("checked"), true).draw();
            });
            /* Event handler */
            $(this.table_element).delegate('tr td:first-child', 'click', function () {
                // Reaction in 1st column was clicked
                let tr = $(this).closest("tr");
                if (tr.hasClass("group")) {
                    return;
                }
                let row = self.datatable.row(tr);
                self.app.dialog_reaction.show(row.data());
            });
            $(this.table_element).on("click", ".cyano-metabolite", function (event) {
                // Any Metabolite in 2nd column was clicked
                const met = self.app.model.metabolite.checked_get("id", this.dataset.id);
                self.app.dialog_metabolite.show(met);
            });
            // First col
            /*table_enzymes.delegate('tr td:first-child', 'click', function() {
                var tr = $(this).closest("tr");
                if (tr.hasClass("group")) {
                    return;
                }
                var row = datatable_enzymes.row(tr);
    
                showEditEnzymeDialog(row.data(), false);
            });*/
            // Enabled in 4th col
            $(this.table_element).delegate('tr td:nth-child(4) input', 'change', function () {
                let row = self.datatable.row($(this).closest("tr"));
                let reaction = row.data();
                reaction.enabled = ($(this).is(":checked"));
                app.command_list.push({
                    "type": "reaction",
                    "op": "edit",
                    "id": reaction.id,
                    "object": {
                        "id": reaction.id,
                        "enabled": reaction.enabled
                    }
                });
                self.invalidate(reaction);
            });
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
        invalidate(reaction) {
            for (let met of reaction.getMetabolites(this.app.model)) {
                this.app.metabolite_page.datatable.row(this.app.model.metabolite.checked_index("id", met.id)).invalidate("data");
            }
            this.datatable.row(this.app.model.reaction.checked_index("id", reaction.id)).invalidate();
        }
    }
    exports.Page = Page;
});
