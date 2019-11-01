function downloadSVG() {
	var s = new XMLSerializer().serializeToString($("#visual > svg")[0]);

	var d = "data:image/svg+xml;base64," + window.btoa(s);

	downloadURI(d, "MBA.svg");
}

function downloadURI(uri, name) {
  var link = document.createElement("a");
  link.download = name;
  link.href = uri;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  delete link;
}

function downloadPNG() {
	var svgelem = $("#visual > svg")[0];
	var w = svgelem.getBBox().width;
	var h = svgelem.getBBox().height;
    var can      = document.createElement('canvas');
    var ctx      = can.getContext('2d');
    var loader   = new Image;                        // Not shown on page

	loader.width  = can.width  = w * 4;
	loader.height = can.height = h * 4;
	loader.onload = function(){
	  	ctx.drawImage( loader, 0, 0, loader.width, loader.height );
		downloadURI(can.toDataURL("image/png"), "MBA.png");
	};
	var svgAsXML = (new XMLSerializer).serializeToString( svgelem );
	loader.src = 'data:image/svg+xml,' + encodeURIComponent( svgAsXML );
}
