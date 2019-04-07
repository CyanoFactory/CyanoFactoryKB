import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

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

<div class="visual_graph">
</div>

<div class="visual_fba">
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

export class Page {
    readonly app: app.AppManager;
    readonly datatable_flux: DataTables.Api;
    readonly source_element: HTMLElement;
    readonly table_element_flux: HTMLElement;

    constructor(where: HTMLElement, app: app.AppManager) {
        this.source_element = where;
        where.appendChild(template.content.cloneNode(true));

        this.table_element_flux = <HTMLElement>where.getElementsByClassName("cyano-flux-list")[0]!;
        this.app = app;

        this.datatable_flux = $(this.table_element_flux).DataTable({
            "deferRender": true,
            "displayLength": 25,
            "order": [[ 1, 'desc' ]],
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
    }

    update() {
    }
}
