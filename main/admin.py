from django.contrib import admin
from .models import Material


class MaterialAdmin(admin.ModelAdmin):
    list_display = ['id', 'type', 'd_1', 'd_2', 'd_3', 'angle']
admin.site.register(Material, MaterialAdmin)
