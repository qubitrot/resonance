$(document).ready(function() {
	$('.sub_func').hide();
	$('.toggle a').click(function() {
		$(this).parent().parent().next().children('td').children('table.sub_func').toggle();
		if ($(this).text() == '+') {
			$(this).text('-');
		} else {
			$(this).text('+');
		}
	});
});
