var http_req;

function toggle_export_popup(id) {
    var div_obj = document.getElementById(id);

    if (div_obj.style.display == "none") {
        div_obj.style.display = "";
    }
    else {
        div_obj.style.display = "none";
    }
}

function close_export_popup(id) {
    var div_obj = document.getElementById(id);

    div_obj.style.display = "none";
}

function export_data(id) {
    //
    // Fetch the marker annotation
    //
    var form_obj = document.getElementById(id + "_frm"); 
    var otype;

    for (i = 0; i < form_obj.otype.length; i++)
        if (form_obj.otype[i].checked)
            otype = form_obj.otype[i].value;

    var url = 
        form_obj.url.value + "&" + 
        "email=" + form_obj.email.value + "&" +
        "otype=" + otype;

    // 
    // Prepare and send XMLHttpRequest Object.
    //
    http_req = false;

    try {
	http_req = new XMLHttpRequest();
    } catch(e) {
	http_req = false;
    }

    if (http_req) {
	http_req.onreadystatechange = process_export;
	http_req.open("GET", url, true);
	http_req.send("");
    }

    var txt_obj = document.getElementById(id + "_txt"); 

    while (txt_obj.childNodes.length > 0)
        txt_obj.removeChild(txt_obj.firstChild); 
}

function gen_export_result(loci, email) {

    var div_obj = document.getElementById("export_popup_txt"); 
    var res_p   = document.createElement("p");

    res_p.innerHTML = "Exporting <span class=\"r\">" + loci + "</span> loci. " +
                      "E-mail will be sent to <span class=\"r\">" + email + "</span>" +
                      " when it is complete.";

    var app_obj = div_obj.appendChild(res_p);

    var res_p       = document.createElement("p");
    res_p.innerHTML = "<a onclick=\"close_export_popup('export_popup')\">close</a>";

    app_obj.appendChild(res_p);
}

function process_export() {
    //
    // Possible readyState values:
    // 0 = uninitialized
    // 1 = loading
    // 2 = loaded
    // 3 = interactive
    // 4 = complete
    //
    if (http_req.readyState == 4) {

        // Check that the status is "OK"
        if (http_req.status == 200) {

            var xml_doc = http_req.responseXML;

            var obj  = xml_doc.getElementsByTagName("loci");
            var loci = obj[0].childNodes[0].nodeValue;

            obj       = xml_doc.getElementsByTagName("email");
            var email = obj[0].childNodes[0].nodeValue;

            obj       = xml_doc.getElementsByTagName("msg");
            var msg   = obj[0].childNodes[0].nodeValue;

            //alert(msg);

            gen_export_result(loci, email);

        } else {
            alert("There was a problem retrieving the XML data:\n" + http_req.statusText);
        }
    }
}