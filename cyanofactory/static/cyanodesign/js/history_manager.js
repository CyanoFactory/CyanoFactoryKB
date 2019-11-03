define(["require", "exports", "./metabolic_model"], function (require, exports, mm) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class HistoryEntry {
        constructor() {
            this.undo = true;
        }
    }
    exports.HistoryEntry = HistoryEntry;
    class HistoryGroup {
    }
    exports.HistoryGroup = HistoryGroup;
    class HistoryManager {
        constructor(history, app) {
            this.history = history;
            this.app = app;
        }
        push(entry, group = null) {
            let hentry = new HistoryEntry();
            hentry.type = entry["type"];
            hentry.op = entry["op"];
            hentry.id = entry["id"];
            hentry.object = entry["object"];
            hentry.group = group;
            let i = 0;
            while (i < this.history.length) {
                if (!this.history[i].undo) {
                    this.history.splice(i, 1);
                }
                else {
                    ++i;
                }
            }
            this.history.push(hentry);
            this.app.history_page.refresh();
        }
        op(idx) {
            const op = this.history[idx]["op"];
            const first = op[0].toUpperCase();
            return first + op.substr(1) + " " + this.history[idx]["type"];
        }
        source(idx) {
            return this.history[idx]["id"];
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
                let s = "";
                let first = true;
                for (const item in obj) {
                    if (!first) {
                        s += ", ";
                    }
                    first = false;
                    s += item + ": " + obj[item];
                }
                return s;
                /*return "Edit reaction " + name;
                if (name != obj["name"]) {
                    return "Renamed reaction " + name + " to " + obj["name"];
                }*/
            }
            else if (op == "delete") {
                return `${id} (${obj["name"]})`;
            }
        }
        apply(idx) {
            const entry = this.history[idx];
            const op = entry.op;
            if (op != "add" && op != "edit" && op != "delete") {
                throw new Error(entry.op + " is not a valid operation!");
            }
            const type = entry.type;
            if (type != "reaction" &&
                type != "metabolite" &&
                type != "compartment" &&
                type != "objective") {
                throw new Error(entry.type + " is not a valid target type!");
            }
            const fn = this.app.model[entry.type];
            const id = entry.id;
            if (op == "add") {
                if (fn.has("id", id)) {
                    throw new Error(entry.type + " " + entry.id + " is already in the model!");
                }
                let obj = fn.create();
                const eobj = entry.object;
                if (id != eobj.id) {
                    throw new Error(entry.type + " " + entry.id + " ID mismatch!");
                }
                for (const key in eobj) {
                    if (eobj.hasOwnProperty(key)) {
                        if (obj.hasOwnProperty(key)) {
                            obj[key] = eobj[key];
                        }
                    }
                }
                fn.add(obj);
            }
            else if (op == "edit") {
                let obj = fn.get("id", id);
                if (obj == null) {
                    throw new Error(entry.type + " " + entry.id + " is not in the model!");
                }
                const eobj = entry.object;
                if (id != eobj.id && fn.has("id", eobj.id)) {
                    throw new Error(entry.type + " " + eobj.id + " is already in the model!");
                }
                for (const key in eobj) {
                    if (eobj.hasOwnProperty(key)) {
                        if (obj.hasOwnProperty(key)) {
                            obj[key] = eobj[key];
                        }
                    }
                }
            }
            else if (op == "delete") {
                if (!fn.has("id", id)) {
                    throw new Error(entry.type + " " + entry.id + " is not in the model!");
                }
                fn.remove("id", id);
            }
        }
        undo(idx) {
            for (let i = idx; i < this.history.length; ++i) {
                this.history[i].undo = false;
            }
            this.app.model = new mm.Model();
            this.app.model.fromJson(this.app.old_models[0]);
            for (let i = 0; i < this.history.length; ++i) {
                if (!this.history[i].undo) {
                    break;
                }
                this.apply(i);
            }
            this.refresh();
        }
        redo(idx) {
            for (let i = 0; i <= idx && i < this.history.length; ++i) {
                if (!this.history[i].undo) {
                    this.history[i].undo = true;
                    this.apply(i);
                }
            }
            this.refresh();
        }
        refresh() {
            this.app.reaction_page.init();
            this.app.metabolite_page.init();
            this.app.settings_page.init();
            this.app.compartment_page.init();
            this.app.history_page.refresh();
        }
    }
    exports.HistoryManager = HistoryManager;
});
