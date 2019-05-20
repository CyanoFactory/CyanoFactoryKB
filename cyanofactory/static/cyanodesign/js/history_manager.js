define(["require", "exports"], function (require, exports) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class HistoryManager {
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
        constructor(history, app) {
            this.history = history;
        }
        op(idx) {
            const op = this.history[idx]["op"];
            const first = op[0].toUpperCase();
            return first + op.substr(1) + " " + this.history[idx]["type"];
        }
        parse(idx) {
            const op = this.history[idx]["op"];
            const typ = this.history[idx]["type"];
            const id = this.history[idx]["id"];
            const obj = this.history[idx]["object"];
            if (op == "add") {
                return `${id} (${obj["name"]})`;
            }
            else if (op == "edit") {
                return "Edit reaction " + name;
                if (name != obj["name"]) {
                    change_list.push("Renamed reaction " + name + " to " + obj["name"]);
                }
            }
            else if (op == "delete") {
                return `${id} (${obj["name"]})`;
            }
        }
        apply(idx) {
        }
        undo(idx) {
        }
    }
    exports.HistoryManager = HistoryManager;
});
