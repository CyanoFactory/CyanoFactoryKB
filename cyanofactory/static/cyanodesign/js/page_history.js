define(["require", "exports", "jquery", "./history_manager", "datatables.net"], function (require, exports, $, history_manager_1) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    let template = document.createElement('template');
    template.innerHTML = `
<table class="cyano-list table table-striped table-hover">
    <thead>
        <tr>
            <th>ID</th>
            <th>Operation</th>
            <th>Target</th>
            <th>Summary</th>
            <th>Undo/Redo</th>
        </tr>
    </thead>
</table>
`;
    class Page {
        constructor(where, app) {
            this.source_element = where;
            where.appendChild(template.content.cloneNode(true));
            this.table_element = where.getElementsByClassName("cyano-list")[0];
            this.app = app;
            const self = this;
            const groupColumn = 5;
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
                        //"width": "5%",
                        "orderable": false,
                        "render": function (data, type, row, meta) {
                            return meta.row + 1;
                        }
                    },
                    {
                        "targets": 1,
                        //"width": "20%",
                        "orderable": false,
                        "render": function (data, type, row, meta) {
                            return self.app.history_manager.op(meta.row);
                        }
                    },
                    {
                        "targets": 2,
                        //"width": "20%",
                        "orderable": false,
                        "render": function (data, type, row, meta) {
                            return self.app.history_manager.source(meta.row);
                        }
                    },
                    {
                        "targets": 3,
                        //"width": "67%",
                        "orderable": false,
                        "render": function (data, type, row, meta) {
                            return self.app.history_manager.parse(meta.row);
                        }
                    },
                    {
                        "targets": 4,
                        //"width": "8%",
                        "orderable": false,
                        "searchable": false,
                        "data": function (data, type, set, meta) {
                            return data;
                        },
                        "render": function (data, type, row, meta) {
                            if (row.undo) {
                                return "<a class='btn btn-default btn-xs undo-button'>Undo</a>";
                            }
                            else {
                                return "<a class='btn btn-default btn-xs redo-button'>Redo</a>";
                            }
                        }
                    },
                    {
                        "targets": groupColumn,
                        "visible": false,
                        "data": function (data, type, set, meta) {
                            return data["group"];
                        }
                    }
                ],
                "order": [[0, "desc"]],
                "displayLength": 100,
                dom: "<'row'<'col-sm-6'l><'col-sm-6'f>>" +
                    "<'row'>" +
                    "<'row'<'col-sm-12'tr>>" +
                    "<'row'<'col-sm-12'B>>" +
                    "<'row'<'col-sm-6'i><'col-sm-6'p>>",
                drawCallback: function (settings) {
                    const api = this.api();
                    const rows = api.rows({ page: 'current' }).nodes();
                    let last = null;
                    api.column(groupColumn, { page: 'current' }).data().each(function (group, i) {
                        const id = group == undefined ? 0 : group["id"];
                        const msg = group == undefined ? "Unsaved changes" : (group["date"] + ": " + group["summary"]);
                        if (last !== id) {
                            $(rows).eq(i).before('<tr class="group"><td colspan="5">' + msg + '</td></tr>');
                            last = id;
                        }
                    });
                }
            });
            // 5th col: undo button clicked
            $(this.table_element).on("click", ".undo-button", function () {
                let row = self.datatable.row($(this).closest("tr"));
                self.app.history_manager.undo(row.index());
            });
            // 5th col: redo button clicked
            $(this.table_element).on("click", ".redo-button", function () {
                let row = self.datatable.row($(this).closest("tr"));
                self.app.history_manager.redo(row.index());
            });
        }
        init(revisions) {
            this.datatable.clear();
            this.datatable.rows.add(this.app.history_manager.history);
            this.datatable.draw();
            if (revisions == null) {
                return;
            }
            let i = 1;
            for (const group of revisions) {
                const hgroup = new history_manager_1.HistoryGroup();
                hgroup.id = i;
                hgroup.date = group["date"];
                hgroup.summary = group["reason"];
                for (const revs of group["changes"]) {
                    this.app.history_manager.push(revs, hgroup);
                }
                ++i;
            }
        }
        refresh() {
            this.init(null);
        }
    }
    exports.Page = Page;
});
