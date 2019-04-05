/*
Copyright (c) 2019 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
*/
var ElementBase = /** @class */ (function () {
    function ElementBase() {
    }
    ElementBase.prototype.greet = function () {
        return "Hello, " + this.name;
    };
    return ElementBase;
}());
