import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";

export class HistoryManager {
    readonly history: any[];

    /*
    app.command_list.push({
        "type": "reaction",
        "op": "edit",
        "id": reaction.id,
        "object": {
            "id": reaction.id,
            "enabled": reaction.enabled
        }
    });
    */

    constructor(history: any[], app: app.AppManager) {
        this.history = history;
    }

    op(idx: number): string {
        const op: string = this.history[idx]["op"];
        const first: string = op[0].toUpperCase();
        return first + op.substr(1) + " " + this.history[idx]["type"];
    }

    parse(idx: number): string {
        const op: string = this.history[idx]["op"];
        const typ: string = this.history[idx]["type"];
        const id: string = this.history[idx]["id"];
        const obj: any = this.history[idx]["object"];

        if (op == "add") {
            return `${id} (${obj["name"]})`;
        } else if (op == "edit") {
            return "Edit reaction " + name;
            if (name != obj["name"]) {
                change_list.push("Renamed reaction " + name + " to " + obj["name"]);
            }
        } else if (op == "delete") {
            return `${id} (${obj["name"]})`;
        }
    }

    apply(idx: number): void {

    }

    undo(idx: number): void {

    }
}
