define(["require", "exports", "jquery", "./history_manager", "datatables.net"], function (require, exports, $, history_manager_1) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    let template = document.createElement('template');
    template.innerHTML = `
<table class="cyano-list table table-striped table-hover">
    <thead>
        <tr>
            <th>Operation</th>
            <th>Summary</th>
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
            this.history_manager = new history_manager_1.HistoryManager(app.command_list);
            const self = this;
            this.datatable = $(this.table_element).DataTable({
                "deferRender": true,
                columns: [
                    {},
                    {}
                ],
                columnDefs: [
                    {
                        "targets": 0,
                        "orderable": false,
                        "render": function (data, type, row, meta) {
                            return self.history_manager.op(meta.row);
                        }
                    },
                    {
                        "targets": 1,
                        "orderable": false,
                        "render": function (data, type, row, meta) {
                            return self.history_manager.parse(meta.row);
                        }
                    }
                ],
                "displayLength": 100,
                dom: "<'row'<'col-sm-6'l><'col-sm-6'f>>" +
                    "<'row'>" +
                    "<'row'<'col-sm-12'tr>>" +
                    "<'row'<'col-sm-12'B>>" +
                    "<'row'<'col-sm-6'i><'col-sm-6'p>>"
            });
        }
        init() {
            this.datatable.clear();
            this.datatable.rows.add(this.app.command_list);
            this.datatable.draw();
        }
        refresh() {
            this.init();
        }
    }
    exports.Page = Page;
});
