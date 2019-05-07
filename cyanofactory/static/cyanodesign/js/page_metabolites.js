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
    let template_filter = document.createElement('template');
    template_filter.innerHTML = `
<div class="col-sm-6">
<label for="cyano-list-filter">Filter metabolites</label>
<select class="cyano-list-filter form-control combobox" multiple="multiple">
    <option selected="selected">Is internal</option>
    <option selected="selected">Is external</option>
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
        //readonly filter_element: HTMLElement;
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
                    else if (options.length === 2) {
                        return 'No filter applied';
                    }
                    else {
                        let labels = [];
                        options.each(function () {
                            labels.push($(this).html());
                        });
                        return labels.join(' ') + '';
                    }
                }
            });
            $.fn.dataTable.ext.search.push(function (settings, data, dataIndex) {
                if (settings.nTable == self.table_element) {
                    const arr = $(where).find(".cyano-list-filter").find("option").map(function () {
                        return this.selected;
                    }).get();
                    const d = self.datatable.data()[dataIndex];
                    const f = [
                        function (e) { return !e.isExternal(); },
                        function (e) { return e.isExternal(); }
                    ];
                    if (!(arr[0] || arr[1])) {
                        return false;
                    }
                    return ((arr[0] ? f[0](d) : !f[0](d)) || (arr[1] ? f[1](d) : !f[1](d)));
                }
                return true;
            });
            where.getElementsByClassName("cyano-regex")[0].addEventListener("click", function () {
                self.datatable.search(self.datatable.search(), $(this).prop("checked"), true).draw();
            });
            /* Event handler */
            // Tooltip on reaction hover
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
            // 1st column clicked
            $(this.table_element).delegate('tr td:first-child', 'click', function () {
                let row = self.datatable.row($(this).closest("tr"));
                app.dialog_metabolite.show(row.data());
            });
            // 4th column: External checkbox
            $(this.table_element).delegate('tr td:nth-child(4) input', 'change', function () {
                let row = self.datatable.row($(this).closest("tr"));
                let metabolite = row.data();
                metabolite.compartment = $(this).is(":checked") ? "e" : "c";
                app.command_list.push({
                    "type": "metabolite",
                    "op": "edit",
                    "id": metabolite.id,
                    "object": {
                        "id": metabolite.id,
                        "name": metabolite.name,
                        "external": metabolite.isExternal()
                    }
                });
                self.invalidate(metabolite);
            });
            /*
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
        /*
            Enzyme.prototype.invalidate = function() {
            this.getMetabolites().forEach(function (obj) {
                datatable_metabolites.row(Metabolite.metabolites.indexOf(obj)).invalidate();
            });
            datatable_enzymes.row(Enzyme.enzymes.indexOf(this)).invalidate();
        };
        Enzyme.prototype.removeFromList = function() {
            var that = this;
            datatable_enzymes.row(function(idx, data, node) {
                return data.name == that.name;
            }).remove();
        };
        Metabolite.prototype.invalidate = function() {
            this.getEnzymes().forEach(function (obj) {
                datatable_enzymes.row(Enzyme.enzymes.indexOf(obj)).invalidate();
            });
            datatable_metabolites.row(Metabolite.metabolites.indexOf(this)).invalidate();
        };
         */
        invalidate(metabolite) {
            for (let reac of metabolite.getReactions()) {
                this.app.reaction_page.datatable.row(this.app.model.reaction.checked_index("id", reac.id)).invalidate("data");
            }
            this.datatable.row(this.app.model.metabolite.checked_index("id", metabolite.id)).invalidate("data");
        }
    }
    exports.Page = Page;
});
