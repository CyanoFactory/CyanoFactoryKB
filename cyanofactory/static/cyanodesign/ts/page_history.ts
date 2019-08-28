import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";
import {HistoryManager} from "./history_manager";

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

export class Page {
    readonly app: app.AppManager;
    readonly datatable: DataTables.Api;
    readonly source_element: HTMLElement;
    readonly table_element: HTMLElement;

    constructor(where: HTMLElement, app: app.AppManager) {
        this.source_element = where;
        where.appendChild(template.content.cloneNode(true));

        this.table_element = <HTMLElement>where.getElementsByClassName("cyano-list")[0]!;
        this.app = app;

        const self: Page = this;

        this.datatable = $(this.table_element).DataTable(<any>{
            "deferRender": true,
            columns: [
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
                    "render": function(data, type, row, meta) {
                        return meta.row + 1;
                    }
                },
                {
                    "targets": 1,
                    //"width": "20%",
                    "orderable": false,
                    "render": function(data, type, row, meta ) {
                        return self.app.history_manager.op(meta.row);
                    }
                },
                {
                    "targets": 2,
                    //"width": "20%",
                    "orderable": false,
                    "render": function(data, type, row, meta ) {
                        return self.app.history_manager.source(meta.row);
                    }
                },
                {
                    "targets": 3,
                    //"width": "67%",
                    "orderable": false,
                    "render": function(data, type, row, meta ) {
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
                    "render": function(data, type, row, meta) {
                        if (row.undo) {
                            return "<a class='btn btn-default btn-xs undo-button'>Undo</a>";
                        } else {
                            return "<a class='btn btn-default btn-xs redo-button'>Redo</a>";
                        }
                    }
                },
            ],
            "order": [[ 0, "desc" ]],
            "displayLength": 100,
            dom:
                "<'row'<'col-sm-6'l><'col-sm-6'f>>" +
                "<'row'>" +
                "<'row'<'col-sm-12'tr>>" +
                "<'row'<'col-sm-12'B>>" +
                "<'row'<'col-sm-6'i><'col-sm-6'p>>"
        });

        // 4th col: undo button clicked
        $(this.table_element).on("click", ".undo-button", function() {
            let row = self.datatable.row($(this).closest("tr"));
            self.app.history_manager.undo(row.index());
        });

        // 4th col: redo button clicked
        $(this.table_element).on("click", ".redo-button", function() {
            let row = self.datatable.row($(this).closest("tr"));
            self.app.history_manager.redo(row.index());
        });
    }

    init() {
        this.datatable.clear();
        this.datatable.rows.add(this.app.history_manager.history);
        this.datatable.draw();
    }

    refresh() {
        this.init();
    }
}
