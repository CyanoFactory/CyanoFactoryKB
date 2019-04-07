define(["require", "exports", "jquery", "datatables.net"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
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
    class Page {
        constructor(where, app) {
            this.source_element = where;
            where.appendChild(template.content.cloneNode(true));
            this.table_element = where.getElementsByClassName("cyano-metabolite-list")[0];
            this.app = app;
            this.datatable = $(this.table_element).DataTable({
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
                        "render": function (data, type, row, meta) {
                            return data.map(function (arg) {
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
                        "render": function (data, type, row, meta) {
                            return data.map(function (arg) {
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
                        render: function (data, type, row, meta) {
                            if (type === 'copy') {
                                return data === true ? "External" : "Internal";
                            }
                            return "<div class='checkbox'> \
                        <input type='checkbox' id='external" + meta.row + "' " + (data ? "checked='checked'" : "") + "> \
                        <label for='external" + meta.row + "'>External</label> \
                        </div>";
                        }
                    }
                ],
                "displayLength": 25,
                dom: "<'row'<'col-sm-6'l><'col-sm-6'f>>" +
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
            // Event handler
            const self = this;
            $(this.table_element).on({
                mouseenter: function () {
                    $(this).data("toggle", "tooltip");
                    $(this).data("placement", "top");
                    let reaction = self.app.model.reaction.checked_get("name", $(this).text());
                    let text = reaction.toString();
                    $(this).prop("title", text);
                    $(this)["tooltip"]("show");
                },
                mouseleave: function () {
                }
            }, ".cyano-enzyme");
        }
        update() {
            this.datatable.clear();
            this.datatable.rows.add(this.app.model.metabolites);
            this.datatable.draw();
        }
    }
    exports.Page = Page;
});
